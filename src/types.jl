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

    LinearRecSystem(arg::T) where {T} = new{T}([], arg, [], [])
end

order(lrs::LinearRecSystem) = length(lrs.mat) - 1
nrows(lrs::LinearRecSystem) = length(lrs.mat) == 0 ? 0 : size(lrs.mat[1], 1)
nfuncs(lrs::LinearRecSystem) = length(lrs.funcs)

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
    # increas order if necessary - add new matrix
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
    push!(lrs.inhom, entry.inhom)
end

lpar(h::Int, d = "") = h == 1 ? "[$(d)" : join(["╭$(d)"; fill("│$(d)", h-2); "╰$(d)"], "\n")
rpar(h::Int, d = "") = h == 1 ? "$(d)]" : join(["$(d)╮"; fill("$(d)│", h-2); "$(d)╯"], "\n")
space(h::Int) = join(fill(" ", h), "\n")
function symstr(h::Int, symbol::String)
    a = fill("   ", h)
    a[Int(ceil(h/2))] = " $(symbol) "
    return join(a, "\n")
end
funcstr(funcs::Vector{String}, arg::String) = join(string.(" ", funcs, "($(arg))"), "\n")

function Base.show(io::IO, lrs::LinearRecSystem)
    h = nrows(lrs)
    lp = lpar(h)
    rp = rpar(h, " ")
    pl = symstr(h, "+")
    eq = symstr(h, "=")

    funcs = string.(lrs.funcs)
    arg = string(lrs.arg)
    m = sprint.(Base.print_matrix, lrs.mat)
    inhom = sprint(Base.print_matrix, lrs.inhom)
    arr = [space(h), lp, m[1], rp, lp, funcstr(funcs, "$(arg)")]
    for i in 2:length(m)
        push!(arr, rp, pl, lp, m[i], rp, lp, funcstr(funcs, "$(arg)+$(i-1)"))
    end
    push!(arr, rp, eq, lp, inhom, rp)

    println(io, "$(typeof(lrs)) of order $(order(lrs)):")
    print(io, mergestr(arr...))
end

function mergestr(strings::String...)
    splits = split.(strings, "\n")
    join(join.(zip(splits...)), "\n")
end