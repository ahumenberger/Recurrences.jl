abstract type Recurrence end
abstract type LinearRecurrence <: Recurrence end

struct CFiniteRecurrence{T} <: LinearRecurrence
    func::T
    arg::T
    coeffs::Vector{T}
    inhom::T

    function CFiniteRecurrence(func::T, arg::T, coeffs::Vector{T}, inhom::T = T(0)) where {T}
        # if any(has.(coeffs, arg))
        #     error("Not a C-finite recurrence.")
        # end
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

coeffs(r::HyperRecurrence) = r.coeffs
order(r::HyperRecurrence) = length(r.coeffs) - 1

function Base.show(io::IO, r::Union{CFiniteRecurrence,HyperRecurrence})
    res = join(["$(c) * $(r.func)($(r.arg + (i - 1)))" for (i, c) in enumerate(r.coeffs)], " + ")
    print(io, "$(res) = $(r.inhom)")
end
