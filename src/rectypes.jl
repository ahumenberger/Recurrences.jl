abstract type Recurrence end
abstract type LinearRecurrence <: Recurrence end

struct HyperRecurrence{T} <: LinearRecurrence
    func::Symbol
    arg::PolyElem{T}
    coeffs::Vector{<:FracElem{<:PolyElem{T}}}
end

# ------------------------------------------------------------------------------

coeffs(r::HyperRecurrence) = r.coeffs
order(r::HyperRecurrence) = length(r.coeffs) - 1

function Base.show(io::IO, r::HyperRecurrence)
    res = join(["$(c) * $(r.func)($(r.arg + (i - 1)))" for (i, c) in enumerate(r.coeffs)], " + ")
    print(io, "$(res) = 0")
end

# ------------------------------------------------------------------------------

function closedform(r::HyperRecurrence{T}; init=nothing) where {T}
    sz = order(r)
    vars = [initvar(string(r.func), i) for i in 0:sz - 1] 
    R, b = PolynomialRing(base_ring(r.arg), vars)
    F = FractionField(R)

    sequences = petkovsek(r.coeffs, r.arg)
    @debug "After Petkovsek" sequences init

    A = [F(s(i)) for i in 0:sz-1, s in sequences]
    b = [F(x) for x in b]
    mA = MatrixSpace(F, size(A)...)(A)
    mB = MatrixSpace(F, length(b), 1)(b)
    mX = Nemo.solve(mA, mB)

    if isnothing(init)
        cfs = Array(mX)
    else
        ini = Array(init)
        cfs = []
        for x in Array(mX)
            n, d = numerator(x), denominator(x)
            n = isone(n) ? one(first(ini)) : sum(i*c for (i, c) in zip(ini, Nemo.coeffs(n)))
            d = isone(d) ? one(first(ini)) : sum(i*c for (i, c) in zip(ini, Nemo.coeffs(d)))
            push!(cfs, n//d)
        end
    end

    F = isnothing(init) ? F : base_ring(base_ring(parent(first(cfs))))
    sequences = [change_coeff_field(F, s) for s in sequences]
    S = parent(first(sequences))
    @debug "Coefficient field" F S

    cform = sum(S(c)*s for (c, s) in zip(cfs, sequences))

    @debug "Closed form" sequences mA mB mX cform
    return cform
end