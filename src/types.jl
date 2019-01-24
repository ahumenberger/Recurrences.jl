export Recurrence, LinearRecurrence, CFiniteRecurrence, LinearRecSystem

abstract type Recurrence end
abstract type LinearRecurrence <: Recurrence end

struct CFiniteRecurrence{T} <: LinearRecurrence
    func::T
    arg::T
    coeffs::Vector{T}
    inhom::T

    function CFiniteRecurrence{T}(func::T, arg::T, coeffs::Vector{T}, inhom::T = T(0)) where {T}
        if !all(is_constant.(coeffs))
            error("Not a C-finite recurrence.")
        end
        new{T}(func, arg, coeffs, inhom)
    end
end

# ------------------------------------------------------------------------------

struct LinearRecSystem{T}
    funcs::Vector{T}
    arg::T
    mat::Vector{Matrix{T}}
    inhom::Vector{T}
end

LinearRecSystem(arg::T) where {T} = LinearRecSystem{T}([], arg, [], [])
LinearRecSystem(funcs::Vector{T}, arg::T, mat::Vector{Matrix{T}}) where {T} = LinearRecSystem{T}(funcs, arg, mat, zeros(T, size(mat[1], 1)))

order(lrs::LinearRecSystem) = length(lrs.mat) - 1
nrows(lrs::LinearRecSystem) = length(lrs.mat) == 0 ? 0 : size(lrs.mat[1], 1)
nfuncs(lrs::LinearRecSystem) = length(lrs.funcs)
ishomogeneous(lrs::LinearRecSystem) = iszero(lrs.inhom)
isdecoupled(lrs::LinearRecSystem) = all(iszero.([m - Diagonal(m) for m in lrs.mat]))

const LinearEntry{T} = NamedTuple{(:coeffs, :inhom), Tuple{Vector{Dict{T,T}}, T}}

findelem(A, e) = findfirst(x -> x == e, A)

function Base.push!(lrs::LinearRecSystem{T}, entry::LinearEntry{T}) where {T}
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

function homogenize!(lrs::LinearRecSystem{T}) where T
    @assert order(lrs) == 1 "Not a recurrence system of order 1."
    if ishomogeneous(lrs)
        return lrs
    end
    push!(lrs.funcs, unique(T))
    n = length(lrs.funcs)

    lrs.mat[1] = vcat(hcat(lrs.mat[1], lrs.inhom * (-1)), zeros(T, 1, n))
    lrs.mat[1][n,n] = -1
    for i in 2:length(lrs.mat)
        lrs.mat[i] = vcat(hcat(lrs.mat[i], zeros(T, n-1)), zeros(T, 1, n))
        lrs.mat[i][n,n] = 1
    end
    
    push!(lrs.inhom, 0)
    fill!(lrs.inhom, 0)
    lrs
end

function decouple(lrs::LinearRecSystem{T}) where T
    if isdecoupled(lrs)
        return lrs
    end
    @assert order(lrs) == 1 "Not a recurrence system of order 1."
    homogenize!(lrs)
    minv = inv(lrs.mat[2])
    matr = [m * minv for m in lrs.mat]
    σ = x -> x |> subs(lrs.arg, lrs.arg+1)
    σinv = x -> x |> subs(lrs.arg, lrs.arg+1)
    δ = x -> σ(x) - x
    C, A = rational_form(copy(-matr[1]), σ, σinv, δ)
    A = -A
    # @info "Zürcher" C A inv(A)

    newlrs = LinearRecSystem(lrs.arg)

    auxcoeff = [C[end,:]; -1] |> subs(lrs.arg, lrs.arg + size(C, 1))
    p = Polynomials.Poly(auxcoeff)

    for i in 1:size(C, 1)
        q = Polynomials.Poly(A[i,:])
        d, r = divrem(p,q)
        @assert r == 0 "Remainder not zero"
        v = Polynomials.coeffs(d)
        k = fill(lrs.funcs[i], length(v))
        d = [Dict([c]) for c in zip(k,v)]
        entry = (coeffs = d, inhom = T(0))
        push!(newlrs, entry)
        # @info "New LRS" newlrs
    end

    newlrs
end

var_count = 0

function unique(::Type{Sym}, n::Int = 1)
    global var_count += n
    if n == 1
        return Sym("ω$var_count")
    end
    return [Sym("ω$i") for i in var_count-n+1:var_count]
end

lpar(h::Int, d = "") = h == 1 ? "[$(d)" : join(["╭$(d)"; fill("│$(d)", h-2); "╰$(d)"], "\n")
rpar(h::Int, d = "") = h == 1 ? "$(d)]" : join(["$(d)╮"; fill("$(d)│", h-2); "$(d)╯"], "\n")
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

function Base.show(io::IO, lrs::LinearRecSystem)
    h = nrows(lrs)
    lp = lpar(h)
    rp = rpar(h, " ")
    pl = symstr(h, "+")
    eq = symstr(h, "=")

    funcs = string.(lrs.funcs)
    arg = string(lrs.arg)
    mstr = sprint.(Base.print_matrix, lrs.mat)
    inhom = sprint(Base.print_matrix, lrs.inhom)
    arr = [space(h), lp, mstr[1], rp, lp, funcstr(funcs, "$(arg)")]

    for i in 2:length(mstr)
        push!(arr, rp, pl, lp, mstr[i], rp, lp, funcstr(funcs, "$(arg)+$(i-1)"))
    end
    push!(arr, rp, eq, lp, inhom, rp)

    println(io, "$(typeof(lrs)) of order $(order(lrs)):")
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