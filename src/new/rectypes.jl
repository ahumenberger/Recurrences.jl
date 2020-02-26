abstract type Recurrence end
abstract type LinearRecurrence <: Recurrence end

struct CFiniteRecurrence{T<:RingElem} <: LinearRecurrence
    func::Symbol
    arg::PolyElem{T}
    coeffs::Vector{T}

    function CFiniteRecurrence(func::Symbol, arg::PolyElem{T}, cfs) where {T<:RingElem}
        ls = T[]
        for c in cfs
            num = numerator(c)
            if isone(denominator(c)) && parent(arg) == parent(num) && length(num) <= 1
                x = base_ring(num)(Nemo.coeff(num, 0))
                push!(ls, x)
            else
                error("Coefficient not an element of the base ring, got $(c)")
            end
        end
        new{T}(func, arg, ls)
    end
end

struct HyperRecurrence{T} <: LinearRecurrence
    func::Symbol
    arg::PolyElem{T}
    coeffs::Vector{<:FracElem{<:PolyElem{T}}}

    function HyperRecurrence(func::Symbol, arg::PolyElem{T}, coeffs) where {T<:RingElem}
        # if any(has.(coeffs, arg))
        #     error("Not a C-finite recurrence.)
        # end
        new{T}(func, arg, coeffs)
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
