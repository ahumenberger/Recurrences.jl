export Recurrence, LinearRecurrence, CFiniteRecurrence, HyperRecurrence

abstract type Recurrence end
abstract type LinearRecurrence <: Recurrence end

struct CFiniteRecurrence{T} <: LinearRecurrence
    func::T
    arg::T
    coeffs::Vector{T}
    inhom::T

    function CFiniteRecurrence(func::T, arg::T, coeffs::Vector{T}, inhom::T = T(0)) where {T}
        if any(has.(coeffs, arg))
            error("Not a C-finite recurrence.")
        end
        new{T}(func, arg, coeffs, inhom)
    end
end

struct HyperRecurrence{T} <: LinearRecurrence
    func::T
    arg::T
    coeffs::Vector{T}
    inhom::T

    function HyperRecurrence(func::T, arg::T, coeffs::Vector{T}, inhom::T = T(0)) where {T}
        # if any(has.(coeffs, arg))
        #     error("Not a C-finite recurrence.")
        # end
        new{T}(func, arg, coeffs, inhom)
    end
end

# ------------------------------------------------------------------------------

coeffs(r::CFiniteRecurrence) = r.coeffs
order(r::CFiniteRecurrence) = length(r.coeffs) - 1

function Base.show(io::IO, r::CFiniteRecurrence)
    res = join(["$(c) * $(r.func)($(r.arg + (i - 1)))" for (i, c) in enumerate(r.coeffs)], " + ")
    print(io, "$(res) = $(r.inhom)")
end



function mroots(poly::Poly{T}) where {T}
    roots = Polynomials.roots(poly)
    Dict([(r, count(x -> x==r, roots)) for r in Base.unique(roots)])
end

function mroots(p::Poly{SymPy.Sym})
    T = SymPy.Sym
    length(p) == 0 && return zeros(T, 0)

    num_leading_zeros = 0
    while p[num_leading_zeros] == zero(T)
        if num_leading_zeros == length(p)-1
            return zeros(T, 0)
        end
        num_leading_zeros += 1
    end
    num_trailing_zeros = 0
    while p[end - num_trailing_zeros] == zero(T)
        num_trailing_zeros += 1
    end
    n = lastindex(p)-(num_leading_zeros + num_trailing_zeros)
    n < 1 && return zeros(T, length(p) - num_trailing_zeros - 1)

    companion = diagm(-1 => ones(T, n-1))
    an = p[end-num_trailing_zeros]
    companion[1,:] = -p[(end-num_trailing_zeros-1):-1:num_leading_zeros] / an

    companion[:eigenvals]()
end