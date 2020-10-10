
mutable struct LinearRecSystem{T<:RingElem}
    funcs::Vector{Symbol}
    arg::T
    mat::Vector{<:MatrixElem{<:FracElem{T}}}
    inhom::MatrixElem{<:FracElem{T}}
end

function LinearRecSystem(arg::T, funcs::Vector{Symbol}=Symbol[]) where {T<:RingElem}
    FF = FractionField(parent(arg))
    sz = length(funcs)
    MS = MatrixSpace(FF, sz, sz)
    LinearRecSystem{T}(funcs, arg, elem_type(MS)[], zero_matrix(FF, 0, 1))
end
# LinearRecSystem(funcs::Vector{S}, arg::S, mat::Vector{Matrix{T}}) where {S,T} = LinearRecSystem{S,T}(funcs, arg, mat, zeros(T, size(mat[1], 1)))

order(lrs::LinearRecSystem) = length(lrs.mat) - 1
nrows(lrs::LinearRecSystem) = length(lrs.mat) == 0 ? 0 : size(lrs.mat[1], 1)
nfuncs(lrs::LinearRecSystem) = length(lrs.funcs)
elem_ring(lrs::LinearRecSystem) = FractionField(parent(lrs.arg))
ishomogeneous(lrs::LinearRecSystem) = iszero(lrs.inhom)
isdecoupled(lrs::LinearRecSystem) = all(i != j ? iszero(m[i,j]) : true for m in lrs.mat for i in 1:nrows(lrs) for j in 1:nrows(lrs))

const LinearEntry{T<:RingElem} = NamedTuple{(:coeffs, :inhom), Tuple{Vector{Dict{Symbol,T}}, T}}

findelem(A, e) = findfirst(x -> x == e, A)

function Base.push!(lrs::LinearRecSystem, entries...)
    R = elem_ring(lrs)
    for entry in entries
        # get all function symbols not occurring in lrs
        newfuncs = setdiff(Iterators.flatten(keys.(entry.coeffs)), lrs.funcs)
        # resize matrices
        newfuncslen = length(newfuncs)
        # @info "" newfuncslen
        sz = nfuncs(lrs) + newfuncslen
        # MS = MatrixSpace(R, sz, sz)
        if newfuncslen > 0
            append!(lrs.funcs, newfuncs)
            # @info "" lrs.funcs newfuncs
            for (i, m) in enumerate(lrs.mat)
                # lrs.mat[i] = hcat(m, zeros(R, size(m, 1), newfuncslen))
                mm = zero_matrix(R, sz, sz)
                for i in 1:Nemo.nrows(m)
                    for j in 1:ncols(m)
                        mm[i,j] = m[i,j]
                    end
                end
                lrs.mat[i] = mm
            end
        end
        # increase order if necessary - add new matrix
        orderdiff = length(entry.coeffs) - length(lrs.mat)
        if orderdiff > 0
            for _ in 1:orderdiff
                # push!(lrs.mat, zeros(R, nrows(lrs), nfuncs(lrs)))
                push!(lrs.mat, zero_matrix(R, sz, sz))
            end
        end
        # add entries to matrices
        rk = maximum(rank(m) for m in lrs.mat)
        for (i, dict) in enumerate(entry.coeffs)
            row = zeros(R, length(lrs.funcs))
            for (key, val) in dict
                # key should be contained in lrs.funcs
                idx = findelem(lrs.funcs, key)
                @assert idx !== nothing
                lrs.mat[i][rk+1, idx] = val
                # row[idx] = val
            end
            # lrs.mat[i] = vcat(lrs.mat[i], transpose(row))
        end
        row = zeros(R, length(lrs.funcs))
        for i in length(entry.coeffs)+1:length(lrs.mat)
            lrs.mat[i] = vcat(lrs.mat[i], transpose(row))
        end
        lrs.inhom = vcat(lrs.inhom, R[entry.inhom;])
    end
    # @info "LinearRecSystem" lrs.funcs lrs.mat lrs.inhom
end


# function firstorder(lrs::LinearRecSystem)
#     order(lrs) == 1 && return lrs

#     R = elem_ring(lrs)
#     n = size(lrs.mat[1], 1)
#     l = order(lrs)
#     funcs = [variables(Symbol, n=(l-1)*n); lrs.funcs]
#     Z = zeros(R, (l-1)*n, n)
#     I = eye(R, (l-1)*n)
#     M = hcat(Z, I)
#     N = hcat(lrs.mat...)
#     @info "" n l funcs Z I M N
#     mat = vcat(M, N)
#     @info "" funcs mat
#     LinearRecSystem{typeof(lrs.arg)}(funcs, lrs.arg, [mat], lrs.inhom)
# end

Base.copy(lrs::LinearRecSystem{T}) where {T} = LinearRecSystem{T}(copy(lrs.funcs), lrs.arg, copy(lrs.mat), copy(lrs.inhom))

function invertible_system!(lrs::LinearRecSystem)
    # We are assuming that lrs is first-order
    M = lrs.mat[1]
    nonzero_cols = [i for i in 1:ncols(M) if !iszero_column(M, i)]
    lrs.mat[1] = M[nonzero_cols, nonzero_cols]
    lrs.mat[2] = lrs.mat[2][nonzero_cols, nonzero_cols]
    lrs.inhom = lrs.inhom[nonzero_cols,:]
    lrs.funcs = lrs.funcs[nonzero_cols]
    # TODO return removed indices
    @debug "Invertible system" nonzero_cols lrs
end

function homogenize(lrs::LinearRecSystem)
    homogenize!(copy(lrs))
end

function homogenize!(lrs::LinearRecSystem)
    @assert order(lrs) == 1 "Not a recurrence system of order 1."
    if ishomogeneous(lrs)
        return lrs
    end

    T = elem_ring(lrs)
    n = length(lrs.funcs) + 1

    lrs.mat[1] = vcat(hcat(lrs.mat[1], lrs.inhom * (-1)), zero_matrix(T, 1, n))
    lrs.mat[1][n,n] = T(-1)
    for i in 2:length(lrs.mat)
        lrs.mat[i] = vcat(hcat(lrs.mat[i], zero_matrix(T, n-1, 1)), zero_matrix(T, 1, n))
        lrs.mat[i][n,n] = T(1)
    end
    
    push!(lrs.funcs, variables(Symbol))
    lrs.inhom = zero_matrix(T, n, 1)
    lrs
end

function monic(lrs::LinearRecSystem)
    m = lrs.mat[end]
    if iszero(m - identity_matrix(elem_ring(lrs), nrows(lrs)))
        return lrs
    end
    minv = inv(m)
    mat = [minv * m for m in lrs.mat]
    inhom = minv * lrs.inhom
    LinearRecSystem(lrs.funcs, lrs.arg, mat, inhom)
end

function blockdiagonal(m::MatrixElem)
    @assert size(m, 1) == size(m, 2) "Matrix is not square."
    R = base_ring(m)
    s = size(m, 1)
    row = zeros(R, 1, s)
    blocks = MatrixElem[]
    i0 = 1
    for i in 1:s-1
        row[1,i] = R(0)
        row[1,i+1] = R(1)
        if Array(m[i,:]) != row
            push!(blocks, m[i0:i, i0:i])
            i0 = i+1
        end
    end
    if i0 <= s
        push!(blocks, m[i0:s, i0:s])        
    end
    # @assert cat(blocks..., dims=(1,2)) == m "Matrix not in block diagonal form."
    blocks
end

(x::FracElem)(y...) = evaluate(numerator(x), y...) // evaluate(denominator(x), y...)
Nemo.AbstractAlgebra.change_base_ring(R::Nemo.AbstractAlgebra.Ring, x::FracElem) = change_base_ring(R, numerator(x)) // change_base_ring(R, denominator(x))

function decouple(lrs::LinearRecSystem)
    @assert order(lrs) == 1 "Not a recurrence system of order 1."
    @assert ishomogeneous(lrs) "Not a homogeneous recurrence system ."
    I = identity_matrix(elem_ring(lrs), nrows(lrs))
    σ = x -> x(lrs.arg-1)
    σinv = x -> x(lrs.arg+1)
    δ = x -> σ(x) - x
    M = map(σ, (-lrs.mat[1])) - I
    C, A = rational_form(copy(M), σ, σinv, δ)
    # C = simplify.(C)
    # A = simplify.(A)
    # @info "" factorize(A)
    @debug "Zürcher" input=-lrs.mat[1] C A M map(σinv, C) (inv(A) * M * map(σ, A) + inv(A) * map(δ, A))
    # @assert inv(A) * M * map(σ, A) + inv(A) * map(δ, A) == C "Zürcher wrong"

    map(σinv, C), A
end

Base.promote_rule(::Type{Frac{T}}, ::Type{Frac{T}}) where T <: RingElem = Frac{T}

function solveblock(C::MatrixElem, initvec::MatrixElem, arg::PolyElem)
    csize = size(C, 1)
    @debug "Solve companion block" C Array(initvec)
    if isone(csize)
        R = base_ring(C)
        coeffs = [-(first(C) + one(R)), one(R)]
    else
        R, x = PolynomialRing(base_ring(C), "x")
        cp1 = sum(sum(c * (-1)^(i-1) * p * x^(j-1) for (j, p) in enumerate(pascal(i-1, alt = true))) for (i, c) in enumerate(Array(C)[end,:]))
        cp2 = sum(p * x^(j-1) for (j, p) in enumerate(pascal(csize, alt = true)))
        coeffpoly = cp2 - cp1
        coeffs = [Nemo.coeff(coeffpoly, i) for i in 0:degree(coeffpoly)]
        @debug "Coeff poly construction for companion block size > 1" degree(coeffpoly) C cp1 cp2
    end
    @debug "Coefficients for recurrence" coeffs 
    # if any(length(numerator(c)) > 1 || length(denominator(c)) > 1 for c in coeffs)
    #     @error "Only C-finite recurrences supported by now"
    #     RecurrenceT = HyperRecurrence
    #     ClosedFormT = HyperClosedForm
    # else
    #     RecurrenceT = CFiniteRecurrence
    #     ClosedFormT = CFiniteClosedForm
    # end
    rec = HyperRecurrence(variables(Symbol), arg, coeffs)
    @debug "Reference recurrence" rec
    cf = closedform(rec, init = initvec)

    # @info "" collect(geometric_sequences(cf))
    # initsubs = Dict(zip(cf.initvec, Array(initvec)[1:length(cf.initvec)]))
    # @debug "Rules for substitution of initial values" initsubs
    # cf = init(cf, initsubs)

    # _arg = change_base_ring(base_ring(base_ring(cf)), arg)
    # @info "" typeof(arg) base_ring(_arg) base_ring(_arg+1) base_ring(cf)
    cforms = Seq[cf]
    for i in 1:csize-1
        push!(cforms, cforms[i](arg+1) - cforms[i](arg))
    end

    cforms
end

const ClosedForm = Pair{Symbol,Seq}

function solve(lrs::LinearRecSystem)
    cforms = Pair{Symbol,Seq}[]
    if isdecoupled(lrs) && ishomogeneous(lrs)
        for i in 1:nrows(lrs)
            coeffs = [m[i,i] for m in lrs.mat]
            rec = HyperRecurrence(lrs.funcs[i], lrs.arg, coeffs)
            cf = closedform(rec)
            push!(cforms, lrs.funcs[i] => cf)
        end
    else
        lrs = monic(lrs)
        @debug "Monic LRS" lrs
        lrs, oldlrs = homogenize(lrs), lrs
        M, A = decouple(lrs)
        blocks = blockdiagonal(M)

        vars = [initvar(string(f), 0) for f in oldlrs.funcs]
        R, initvec = PolynomialRing(base_ring(lrs.arg), vars)
        F = FractionField(R)

        if nfuncs(oldlrs) != nfuncs(lrs)
            push!(initvec, one(R))
        end

        _arg = change_base_ring(F, lrs.arg)
        M = map(x->change_base_ring(F, x), lrs.mat[1])
        S = FractionField(parent(_arg))
        MS = MatrixSpace(S, size(A)...)
        # _mat = map(m->change_base_ring(F, m), lrs.mat)
        _initvec = MatrixSpace(S, length(initvec), 1)([S(x) for x in initvec])
        _A = map(x->change_base_ring(F, x), A)
        # @info "" _arg |> typeof  M

        if lrs.funcs != oldlrs.funcs
            # Assume lrs got homogenized, therefore initial value of introduced variable is 1
            _initvec[length(_initvec), 1] = F(1)
        end
        maxsize = maximum(size.(blocks, 1))
        # M = _mat[1]
        J = identity_matrix(M)
        D = -M
        imat = zero_matrix(S, size(M, 1), 0)
        for i in 0:maxsize-1
            ivec = J * _initvec
            imat = hcat(imat, ivec)
            J = MS(Array(map(x->x(i+1), D))) * J
        end
        @debug "Matrix containing initial values" imat inv(_A)*imat
        result = Seq[]
        i0 = 1
        for C in blocks
            append!(result, solveblock(C, (inv(_A)*imat)[i0,:], lrs.arg))
            i0 += size(C, 1)
        end
        sR = parent(first(result))
        rs = MatrixSpace(sR, length(result), 1)(result)
        result = map(sR, _A) * rs
        @debug "Final closed forms" result

        for (i, v) in enumerate(oldlrs.funcs)
            push!(cforms, v => Array(result)[i])
        end
    end
    cforms
end

# var_count = 0

# function unique(::Type{Sym}, n::Int = 1)
#     global var_count += n
#     if n == 1
#         return Sym("ω$var_count")
#     end
#     return [Sym("ω$i") for i in var_count-n+1:var_count]
# end

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
    mstr = [sprint(Base.print_matrix, Array(m)) for m in lrs.mat]
    inhom = sprint(Base.print_matrix, Array(lrs.inhom))
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