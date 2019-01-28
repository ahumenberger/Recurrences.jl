export Recurrence, LinearRecurrence, CFiniteRecurrence

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