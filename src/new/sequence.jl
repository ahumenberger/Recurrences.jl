mutable struct SeqRing{T <: PolyElem} <: Ring
    base_ring::PolyRing
 
    function SeqRing{T}(R::PolyRing, cached::Bool = true) where T <: PolyElem
       if cached && haskey(SequenceID, R)
          return SequenceID[R]::SeqRing{T}
       else
          z = new{T}(R)
          if cached
            SequenceID[R] = z
          end
          return z
       end
    end
end
 
const SequenceID = Dict{Ring, Ring}()

function SequenceRing(R::Ring, s::AbstractString; cached::Bool = true)
    S, v = PolynomialRing(R, s)
    parent_obj = SeqRing{elem_type(S)}(S, cached)
    return parent_obj, v
end

SequenceRing(R::PolyRing; cached::Bool = true) = SeqRing{elem_type(R)}(R, cached)

# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------


 mutable struct Seq{T <: FieldElem} <: RingElem
    terms::Vector{HyperTerm{T}}
    parent::SeqRing{<:PolyElem{T}}
 
    Seq{T}() where T <: FieldElem = new{T}(Array{HyperTerm{T}}(undef, 0))
 
    function Seq{T}(b::Vector{HyperTerm{T}}) where T <: FieldElem
       z = new{T}(b)
       return z
    end
 
    # Seq{T}(a::T) where T <: PolyElem = iszero(a) ? new{T}(Array{HyperTerm{T}}(undef, 0)) : new{T}([HyperTerm{T}(a)], 1)
 end



# function Base.show(io::IO, t::Seq)
#     if t.length == 0
#         show(io, 0)
#     else
#         show(io, t.terms[1])
#         for term in t.terms[2:end]
#             show(io, " + ", term)
#         end
#     end
# end

# # ------------------------------------------------------------------------------

parent_type(::Type{Seq{T}}) where T <: FieldElem = SeqRing{T}

elem_type(::Type{SeqRing{T}}) where T <: PolyElem = Seq{T}

base_ring(R::SeqRing{T}) where T <: PolyElem = R.base_ring::parent_type(T)
base_ring(a::Seq) = base_ring(parent(a))

parent(a::Seq) = a.parent

isdomain_type(::Type{Seq}) = false
isexact_type(::Type{Seq}) = true

Base.iszero(a::Seq) = iszero(length(a.terms))

Base.one(R::SeqRing) = R(1)
Base.zero(R::SeqRing) = R()

Base.one(s::Seq) = one(parent(s))
Base.zero(s::Seq) = zero(parent(s))

function Base.hash(a::Seq, h::UInt)
    b = 0x53dd43cd511044d1%UInt
    for t in a.terms
        b = xor(b, xor(hash(t, h), h))
        b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
    end
    return b
end

function Base.isequal(x::Seq{T}, y::Seq{T}) where {T <: FieldElem}
    parent(x) != parent(y) && return false
    length(x.terms) != length(y.terms) && return false
    for i in 1:length(x.terms)
       x.terms[i] != y.terms[i] && return false
    end
    return true
end

Base.:(==)(x::Seq{T}, y::Seq{T}) where {T <: FieldElem} = isequal(x, y)

# ------------------------------------------------------------------------------

function (a::SeqRing{T})() where {T <: PolyElem}
    S = elem_type(base_ring(base_ring(a)))
    z = Seq{S}()
    z.parent = a
    return z
end

function (a::SeqRing{T})(b::Union{Integer, Rational}) where {T <: PolyElem}
    F = FractionField(base_ring(a))
    return a(F(b))
end

function (a::SeqRing{T})(b::S) where {T <: PolyElem, S <: FieldElem}
    parent(b) != base_ring(base_ring(a)) && error("Unable to coerce to sequence")
    F = FractionField(base_ring(a))
    return a(F(b))
end

function (a::SeqRing{T})(b::T) where {T <: PolyElem}
    parent(b) != base_ring(a) && error("Unable to coerce to sequence")
    F = FractionField(parent(b))
    return a(F(b))
end

function (a::SeqRing{T})(b::FracElem{T}) where {T <: PolyElem}
    base_ring(b) != base_ring(a) && error("Unable to coerce to sequence")
    S = elem_type(base_ring(base_ring(a)))
    t = HyperTerm{S}(b)
    z = Seq{S}([t])
    z.parent = a
    return z
end

function (a::SeqRing{T})(b::Union{Integer, Rational}, c::T) where {T <: PolyElem}
    F = FractionField(base_ring(a))
    return a(F(b), c)
end

function (a::SeqRing{T})(b::S, c::T) where {T <: PolyElem, S <: FieldElem}
    parent(b) != base_ring(base_ring(a)) && error("Unable to coerce to sequence")
    F = FractionField(base_ring(a))
    return a(F(b), c)
end

function (a::SeqRing{T})(b::T, c::T) where {T <: PolyElem}
    parent(b) != base_ring(a) && error("Unable to coerce to sequence")
    F = FractionField(parent(b))
    return a(F(b), c)
end

function (a::SeqRing{T})(b::FracElem{T}, c::T) where {T <: PolyElem}
    base_ring(b) != base_ring(a) && error("Unable to coerce to sequence")
    S = elem_type(base_ring(base_ring(a)))
    t = HyperTerm{S}(b, c)
    z = Seq{S}([t])
    z.parent = a
    return z
end

# ------------------------------------------------------------------------------

function Base.show(io::IO, a::Seq)
    if iszero(a)
        show(io, 0)
    else
        show(io, a.terms[1])
        for t in a.terms[2:end]
            print(io, " + ")
            show(io, t)
        end
    end
end

# ------------------------------------------------------------------------------

function check_parent(a::Seq, b::Seq, throw::Bool = true)
    c = parent(a) != parent(b)
    c && throw && error("Incompatible polynomial rings in sequence operation")
    return !c
 end

# ------------------------------------------------------------------------------

function Base.:-(f::Seq{T}) where {T <: FieldElem}
    z = Seq{T}([-t for t in f.terms])
    z.parent = f.parent
    return z
end

function Base.:+(f::Seq{T}, g::Seq{T}) where T <: FieldElem
    check_parent(f, g)
    terms = f.terms
    for t in g.terms
        idx = findfirst(x->(issummable(x, t)), terms)
        if idx === nothing
            push!(terms, t)
        else
            terms[idx] += t
        end

    end
    filter!(!iszero, terms)
    isempty(terms) && return zero(f)
    z = Seq{T}(terms)
    z.parent = f.parent
    return z
end

function Base.:-(f::Seq{T}, g::Seq{T}) where T <: FieldElem
    check_parent(f, g)
    terms = f.terms
    for t in g.terms
        idx = findfirst(x->(issummable(x, t)), terms)
        if idx === nothing
            push!(terms, -t)
        else
            terms[idx] -= t
        end
    end
    filter!(!iszero, terms)
    isempty(terms) && return zero(f)
    z = Seq{T}(terms)
    z.parent = f.parent
    return z
end

function Base.:*(f::Seq{T}, g::Seq{T}) where T <: FieldElem
    check_parent(f, g)
    seqs = Seq{T}[]
    for a in f.terms
        for b in g.terms
            z = Seq{T}([a*b])
            z.parent = f.parent
            push!(seqs, z)
        end
    end
    sum(seqs)
end

function Base.:^(f::Seq{T}, i::Integer) where T <: FieldElem
    isnegative(i) && error("Only nonnegative powers allowed")
    iszero(i) && return one(f)
    prod(f for _ in 1:i)
end

# ------------------------------------------------------------------------------

function (a::Seq{T})(b::Integer) where T <: FieldElem
    sum(t(b) for t in a.terms)
end

function (a::Seq{T})(b::PolyElem{T}) where T <: FieldElem
    base_ring(a) != parent(b) && error("Polynomial rings do not match")
    (Nemo.degree(b) != 1 || lead(b) != 1 || denominator(Nemo.coeff(b, 0)) != 1) && error("Argument must be of form n+c where c is an integer")
    z = Seq{T}([t(b) for t in a.terms])
    z.parent = a.parent
    return z
end