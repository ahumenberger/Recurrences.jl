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
    print(io, "$(res) = $(r.inhom)")
end

# ------------------------------------------------------------------------------

function closedform(r::HyperRecurrence{T}) where {T}
    sz = order(r)
    vars = [initvar(string(r.func), i) for i in 0:sz - 1] 
    R, b = PolynomialRing(base_ring(r.arg), vars)
    F = FractionField(R)

    sequences = petkovsek(r.coeffs, r.arg)

    A = [F(s(i)) for i in 0:sz-1, s in sequences]
    b = [F(x) for x in b]
    mA = MatrixSpace(F, size(A)...)(A)
    mB = MatrixSpace(F, length(b), 1)(b)
    mX = Nemo.solve(mA, mB)

    seqs = [change_coeff_field(F, s) for s in sequences]
    S = parent(first(seqs))

    cform = sum(S(c)*s for (c, s) in zip(Array(mX), seqs))

    @debug "Closed form" seqs sequences mA mB mX  cform
    return cform
end