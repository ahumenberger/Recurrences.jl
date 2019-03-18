export LinearRecSystem
export decouple, homogenize!, solve

struct LinearRecSystem{T}
    funcs::Vector{T}
    arg::T
    mat::Vector{Matrix{T}}
    inhom::Vector{T}
end

LinearRecSystem(arg::T, funcs::Vector{T}=T[]) where {T} = LinearRecSystem{T}(funcs, arg, [], [])
LinearRecSystem(funcs::Vector{T}, arg::T, mat::Vector{Matrix{T}}) where {T} = LinearRecSystem{T}(funcs, arg, mat, zeros(T, size(mat[1], 1)))

order(lrs::LinearRecSystem) = length(lrs.mat) - 1
nrows(lrs::LinearRecSystem) = length(lrs.mat) == 0 ? 0 : size(lrs.mat[1], 1)
nfuncs(lrs::LinearRecSystem) = length(lrs.funcs)
ishomogeneous(lrs::LinearRecSystem) = iszero(lrs.inhom)
isdecoupled(lrs::LinearRecSystem) = all(iszero.([m - Diagonal(m) for m in lrs.mat]))

const LinearEntry{T} = NamedTuple{(:coeffs, :inhom), Tuple{Vector{Dict{T,T}}, T}}

findelem(A, e) = findfirst(x -> x == e, A)

function Base.push!(lrs::LinearRecSystem{T}, entries::LinearEntry{T}...) where {T}
    for entry in entries
        # get all function symbols not occurring in lrs
        newfuncs = setdiff(Iterators.flatten(keys.(entry.coeffs)), lrs.funcs)
        # resize matrices
        newfuncslen = length(newfuncs)
        # @info "" newfuncslen
        if newfuncslen > 0
            append!(lrs.funcs, newfuncs)
            # @info "" lrs.funcs newfuncs
            for (i, m) in enumerate(lrs.mat)
                lrs.mat[i] = hcat(m, zeros(T, size(m, 1), newfuncslen))
            end
        end
        # increase order if necessary - add new matrix
        orderdiff = length(entry.coeffs) - length(lrs.mat)
        if orderdiff > 0
            for _ in 1:orderdiff
                push!(lrs.mat, zeros(T, nrows(lrs), nfuncs(lrs)))
            end
        end
        # add entries to matrices
        for (i, dict) in enumerate(entry.coeffs)
            row = zeros(T, length(lrs.funcs))
            for (key, val) in dict
                # key should be contained in lrs.funcs
                idx = findelem(lrs.funcs, key)
                # @info lrs.funcs
                @assert idx != nothing
                row[idx] = val
            end
            # @info "Test" lrs.mat[i] row
            lrs.mat[i] = vcat(lrs.mat[i], transpose(row))
        end
        row = zeros(T, length(lrs.funcs))
        for i in length(entry.coeffs)+1:length(lrs.mat)
            lrs.mat[i] = vcat(lrs.mat[i], transpose(row))
        end
        push!(lrs.inhom, entry.inhom)
    end
    # @info "LinearRecSystem" lrs.funcs lrs.mat lrs.inhom
end

function firstorder(lrs::LinearRecSystem{T}) where {T}
    if order(lrs) == 1
        return lrs
    end

    n = size(lrs.mat[1], 1)
    l = length(lrs.mat)
    funcs = [unique(T, (l-1)*n); lrs.funcs]
    Z = zeros(T, (l-1)*n, n)
    I = eye(T, (l-1)*n)
    M = hcat(Z,I)
    N = hcat(lrs.mat...)
    mat = vcat(M, N)
    LinearRecSystem(funcs, lrs.arg, [mat], lrs.inhom)
end

Base.copy(lrs::LinearRecSystem) = LinearRecSystem(copy(lrs.funcs), lrs.arg, copy(lrs.mat), copy(lrs.inhom))

function homogenize(lrs::LinearRecSystem)
    homogenize!(copy(lrs))
end

function homogenize!(lrs::LinearRecSystem{T}) where {T}
    @assert order(lrs) == 1 "Not a recurrence system of order 1."
    if ishomogeneous(lrs)
        return lrs
    end

    n = length(lrs.funcs) + 1

    lrs.mat[1] = vcat(hcat(lrs.mat[1], lrs.inhom * (-1)), zeros(T, 1, n))
    lrs.mat[1][n,n] = -1
    for i in 2:length(lrs.mat)
        lrs.mat[i] = vcat(hcat(lrs.mat[i], zeros(T, n-1)), zeros(T, 1, n))
        lrs.mat[i][n,n] = 1
    end
    
    push!(lrs.funcs, variables(T))
    push!(lrs.inhom, 0)
    fill!(lrs.inhom, 0)
    lrs
end

function monic(lrs::LinearRecSystem)
    if iszero(lrs.mat[end] - UniformScaling(1))
        return lrs
    end
    minv = inv(lrs.mat[end])
    mat = [minv * m for m in lrs.mat]
    inhom = minv * lrs.inhom
    LinearRecSystem(lrs.funcs, lrs.arg, mat, inhom)
end

function blockdiagonal(m::Matrix{T}) where {T}
    @assert size(m, 1) == size(m, 2) "Matrix is not square."
    s = size(m, 1)
    row = zeros(T, s)
    blocks = Matrix{T}[]
    i0 = 1
    for i in 1:s-1
        row[i] = 0
        row[i+1] = 1
        if m[i,:] != row
            push!(blocks, m[i0:i, i0:i])
            i0 = i+1
        end
    end
    if i0 <= s
        push!(blocks, m[i0:s, i0:s])        
    end
    @assert cat(blocks..., dims=(1,2)) == m "Matrix not in block diagonal form."
    blocks
end

function decouple(lrs::LinearRecSystem{T}) where {T}
    @assert order(lrs) == 1 "Not a recurrence system of order 1."
    @assert ishomogeneous(lrs) "Not a homogeneous recurrence system ."

    σ = x -> subs(x, lrs.arg => lrs.arg-1)
    σinv = x -> subs(x, lrs.arg => lrs.arg+1)
    δ = x -> σ(x) - x
    M = σ.(-lrs.mat[1]) - UniformScaling(1)
    C, A = rational_form(copy(M), σ, σinv, δ)
    C = simplify.(C)
    A = simplify.(A)

    @debug "Zürcher" input=-lrs.mat[1] C A M simplify.(inv(A) * M * A)
    @assert simplify.(inv(A) * M * A) == C "Zürcher wrong"

    σinv.(C), A
end

function solveblock(C::Matrix{T}, initvec::Vector{T}, arg::T) where {T}
    csize = size(C, 1)
    @debug "Solve companion block" C
    coeffpoly = sum(c * Poly(pascal(i-1, alt = true)) for (i, c) in enumerate(C[end,:])) + Poly(pascal(csize, alt = true))
    @debug "Coefficients for recurrence" simplify.(Polynomials.coeffs(coeffpoly))
    coeffs = Polynomials.coeffs(coeffpoly)
    if any(arg in free_symbols(c) for c in coeffs)
        @error "Only C-finite recurrences supported by now"
        RecurrenceT = HyperRecurrence
        ClosedFormT = HyperClosedForm
    else
        RecurrenceT = CFiniteRecurrence
        ClosedFormT = CFiniteClosedForm
    end
    rec = RecurrenceT(variables(T), arg, coeffs)
    @debug "Reference recurrence" rec
    cf = closedform(rec)

    initsubs = Dict(zip(cf.initvec, initvec[1:length(cf.initvec)]))
    @debug "Rules for substitution of initial values" initsubs
    cf = init(cf, initsubs)

    cforms = ClosedForm[cf]
    for i in 1:csize-1
        push!(cforms, cforms[i](arg+1) - cforms[i](arg))
    end

    cforms
end

function solve(lrs::LinearRecSystem{T}) where {T}
    cforms = CFiniteClosedForm[]
    if isdecoupled(lrs) && ishomogeneous(lrs)
        for i in 1:nrows(lrs)
            coeffs = [m[i,i] for m in lrs.mat]
            rec = CFiniteRecurrence(lrs.funcs[i], lrs.arg, coeffs)
            cf = closedform(rec)
            push!(cforms, cf)
        end
    else
        lrs = monic(lrs)
        @debug "Monic LRS" lrs
        lrs, oldlrs = homogenize(lrs), lrs
        M, A = decouple(lrs)
        blocks = blockdiagonal(M)

        initvec = [initvariable(f, 0) for f in lrs.funcs]
        if lrs.funcs != oldlrs.funcs
            # Assume lrs got homogenized, therefore initial value of introduced variable is 1
            initvec[end] = T(1)
        end
        maxsize = maximum(size.(blocks, 1))
        M = -lrs.mat[1]
        S = UniformScaling(1)
        D = inv(A)
        imat = Matrix{T}(undef, size(M, 1), 0)
        for i in 0:maxsize-1
            ivec = (subs(D, lrs.arg => i) * S * initvec)
            imat = hcat(imat, ivec)
            S = subs(M, lrs.arg, i+1) * S
        end
        @debug "Matrix containing initial values" imat
        result = ClosedForm[]
        i0 = 1
        for C in blocks
            append!(result, solveblock(C, imat[i0,:], lrs.arg))
            i0 += size(C, 1)
        end
        result = A * result
        @debug "Final closed forms" result

        for i in 1:length(oldlrs.funcs)
            push!(cforms, ClosedForm(lrs.funcs[i], result[i]))
        end
    end
    cforms
end

initvariable(v::T, i::Union{T, Int64}) where {T} = T("$(string(v))$(i)$(i)")

var_count = 0

function unique(::Type{Sym}, n::Int = 1)
    global var_count += n
    if n == 1
        return Sym("ω$var_count")
    end
    return [Sym("ω$i") for i in var_count-n+1:var_count]
end

# lpar(h::Int, d = "") = h == 1 ? "[$(d)" : join(["╭$(d)"; fill("│$(d)", h-2); "╰$(d)"], "\n")
# rpar(h::Int, d = "") = h == 1 ? "$(d)]" : join(["$(d)╮"; fill("$(d)│", h-2); "$(d)╯"], "\n")
lpar(h::Int, d = "") = h == 1 ? "($(d)" : join(["⎛$(d)"; fill("⎜$(d)", h-2); "⎝$(d)"], "\n")
rpar(h::Int, d = "") = h == 1 ? "$(d))" : join(["$(d)⎞"; fill("$(d)⎟", h-2); "$(d)⎠"], "\n")

space(h::Int) = join(fill(" ", h), "\n")
function symstr(h::Int, symbol::String)
    a = fill("   ", h)
    a[Int(ceil(h/2))] = " $(symbol) "
    return join(a, "\n")
end
function funcstr(funcs::Vector{String}, arg::String)
    a = string.(" ", funcs, "($(arg))")
    maxlen = maximum(length.(a))
    for i in 1:length(a)
        a[i] = string(a[i], repeat(" ", maxlen - length(a[i])))
    end
    join(a, "\n") 
end

Base.summary(io::IO, lrs::LinearRecSystem) = println(io, "$(nrows(lrs))-element $(typeof(lrs)) of order $(order(lrs)):")

function Base.show(io::IO, lrs::LinearRecSystem)
    summary(io, lrs)

    h = nrows(lrs)
    if h == 0
        return
    end

    lp = lpar(h)
    rp = rpar(h, " ")
    pl = symstr(h, "+")
    eq = symstr(h, "=")

    funcs = string.(lrs.funcs)
    arg = string(lrs.arg)
    mstr = sprint.(Base.print_matrix, lrs.mat)
    inhom = sprint(Base.print_matrix, lrs.inhom)
    arr = [space(h), lp, mstr[end], rp, lp, funcstr(funcs, "$(arg)+$(length(mstr) - 1)")]

    for i in length(mstr)-1:-1:1
        push!(arr, rp, pl, lp, mstr[i], rp, lp, funcstr(funcs, i == 1 ? arg : "$(arg)+$(i-1)"))
    end
    push!(arr, rp, eq, lp, inhom, rp)

    if h == 1
        print(io, join(arr))
    else
        print(io, mergestr(arr...))
    end
end

function mergestr(strings::String...)
    splits = split.(strings, "\n")
    rows = length(splits[1])
    cols = length(splits)
    matr = reshape(collect(Iterators.flatten(splits)), rows, cols)
    join([join(matr[i,:]) for i in 1:size(matr, 1)], "\n")
end

function Base.show(io::IO, entry::LinearEntry{T}) where {T}
    arg = "n" # random function argument
    arr = String[]
    for (i, c) in enumerate(entry.coeffs)
        for (k, v) in c
            argstr = i == 1 ? arg : string(arg, "+$(i-1)")
            push!(arr, "$(v)*$(k)($(argstr))")
        end
    end
    str = string(" ", join(arr, " + "), " = $(entry.inhom)")
    println(io, "LinearEntry{$(T)} with default argument $(arg):")
    print(io, str)
end