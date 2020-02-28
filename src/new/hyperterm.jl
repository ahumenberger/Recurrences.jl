struct HyperTerm{T <: FieldElem} 
    coeff::FracElem{<: PolyElem{T}}
    geom::T
    fact::FracElem{<: PolyElem{T}}
    power::PolyElem{T}

    function HyperTerm{T}(a::FracElem{<: PolyElem{T}}, b::PolyElem{T}) where T <: FieldElem
        base_ring(a) != parent(b) && error("Polynomial rings do not match")
        c1, f1 = monic_factors(numerator(a))
        c2, f2 = monic_factors(denominator(a))
        c = Nemo.coeff(c1, 0) // Nemo.coeff(c2, 0)
        p = _prod(f1)//_prod(f2)
        new{T}(one(a), c, parent(a)(p), b)
    end

    function HyperTerm{T}(a::FracElem{<: PolyElem{T}}) where T <: FieldElem
        R = base_ring(a)
        new{T}(a, one(base_ring(R)), one(a), R([0,1]))
    end

    function HyperTerm{T}(a::FracElem{<: PolyElem{T}}, b::T, c::FracElem{<:PolyElem{T}}, d::PolyElem{T}) where T <: FieldElem
        new{T}(a, b, c, d)
    end
end

function normalise(a::HyperTerm{T}) where T <: FieldElem
    norm = parent(a.power)([0, 1])
    a.power == norm && return a

    # TODO: make sure that a.power is of the form n+c where c is an integer
    d = convert(Int, numerator(Nemo.coeff(a.power - norm, 0)))
    c = a.coeff
    if isnegative(d)
        c *= 1//prod(a.fact + i for i in 1:abs(d))
    else
        c *= prod(a.fact - i for i in 0:d-1)
    end
    c *= a.geom^d
    f = isone(a.fact) ? a.fact : a.fact - d

    HyperTerm{T}(c, a.geom, f, norm)
end

Base.isone(a::HyperTerm) = isone(a.coeff) && isone(a.geom) && isone(a.fact)
Base.iszero(a::HyperTerm) = iszero(a.coeff)
Base.one(a::HyperTerm{T}) where T <: FieldElem = HyperTerm{T}(one(a.coeff))
Base.zero(a::HyperTerm{T}) where T <: FieldElem = HyperTerm{T}(zero(a.coeff))

function Base.hash(a::HyperTerm, h::UInt)
    b = 0x53dd43cd511044d1%UInt
    for t in (a.coeff, a.geom, a.fact, a.power)
        b = xor(b, xor(hash(t, h), h))
        b = (b << 1) | (b >> (sizeof(Int)*8 - 1))
    end
    return b
end

function Base.isequal(a::HyperTerm, b::HyperTerm) where T <: FieldElem
    x, y = normalise(a), normalise(b)
    x.coeff == y.coeff && x.geom == y.geom && x.fact == y.fact && x.power == y.power
end

Base.:(==)(x::HyperTerm, y::HyperTerm) = isequal(x, y)

function Base.show(io::IO, a::HyperTerm)
    if iszero(a)
        print(io, 0)
    elseif isone(a)
        print(io, 1)
    else
        out = false
        if !isone(a.coeff)
            print(io, a.coeff)
            out = true
        end
        if !isone(a.geom)
            if out
                print(io, "*")
            end
            q = needs_parentheses(a.geom)
            p = needs_parentheses(a.power)
            print(io, p ? "(" : "", a.geom, p ? ")" : "", "^", p ? "(" : "", a.power, p ? ")" : "")
            out = true
        end
        if !isone(a.fact)
            if out
                print(io, "*")
            end
            p = needs_parentheses(a.fact)
            print(io, p ? "(" : "", a.fact, p ? ")" : "", "^[", a.power, "]")
        end
    end
end

function issummable(f::HyperTerm{T}, g::HyperTerm{T}) where T <: FieldElem
    a, b = normalise(f), normalise(g)
    a.geom == b.geom && a.fact == b.fact
end

Base.:-(a::HyperTerm{T}) where T <: FieldElem = HyperTerm{T}(-a.coeff, deepcopy(a.geom), deepcopy(a.fact), deepcopy(a.power))

function Base.:+(f::HyperTerm{T}, g::HyperTerm{T}) where T <: FieldElem
    a, b = normalise(f), normalise(g)
    !issummable(a, b) && error("Cannod the sum of $f and $g")
    HyperTerm{T}(a.coeff+b.coeff, deepcopy(a.geom), deepcopy(a.fact), deepcopy(a.power))
end

function Base.:-(f::HyperTerm{T}, g::HyperTerm{T}) where T <: FieldElem
    a, b = normalise(f), normalise(g)
    !issummable(a, b) && error("Cannod the sum of $f and $g")
    HyperTerm{T}(a.coeff-b.coeff, deepcopy(a.geom), deepcopy(a.fact), deepcopy(a.power))
end

function Base.:*(f::HyperTerm{T}, g::HyperTerm{T}) where T <: FieldElem
    a, b = normalise(f), normalise(g)
    HyperTerm{T}(a.coeff*b.coeff, a.geom*b.geom, a.fact*b.fact, deepcopy(a.power))
end

function Base.:^(f::HyperTerm{T}, p::Integer) where T <: FieldElem
    a = normalise(f)
    HyperTerm{T}(a.coeff^p, a.geom^p, a.fact^p, deepcopy(a.power))
end

# ------------------------------------------------------------------------------

function (a::HyperTerm{T})(b::Integer) where T <: FieldElem
    z = a.coeff(b)
    z *= a.geom^b
    if !isone(a.fact)
        z *= fallfactorial(a.fact(b), b)
    end
    return z
end

function (a::HyperTerm{T})(b::PolyElem{T}) where T <: FieldElem
    HyperTerm{T}(a.coeff(b), deepcopy(a.geom), a.fact(b), a.power(b))
end