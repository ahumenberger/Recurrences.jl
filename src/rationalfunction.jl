struct RationalFunction{T}
    num::Poly{T}
    den::Poly{T}

    function RationalFunction{T}(num::Poly{T}, den::Poly{T}) where T
        num == den == zero(T) && __throw_rational_argerror(T)
        num2, den2 = divgcd(num, den)
        new(num2, den2)
    end
end
@noinline __throw_rational_argerror(T) = throw(ArgumentError("invalid rational: zero($T)//zero($T)"))

RationalFunction(n::Poly{T}, d::Poly{T}) where {T} = RationalFunction{T}(n, d)
RationalFunction(n::Poly, d::Poly) = RationalFunction(promote(n, d)...)
RationalFunction(n::Poly, d::Number) = RationalFunction(n, Poly(d, n.var))
RationalFunction(n::Number, d::Poly) = RationalFunction(Poly(n, d.var), d)
RationalFunction(n::Poly{T}) where {T} = RationalFunction(n, one(T))
RationalFunction(n::T) where {T} = RationalFunction(n, one(T))
RationalFunction(n::Number, d::Number) = RationalFunction(Poly(n), Poly(d))
RationalFunction(x::RationalFunction) = x

function divgcd(x::Poly, y::Poly)
    g = gcd(x, y)
    div(x, g), div(y, g)
end

Base.://(n::Poly{T}, d::Poly{T}) where T = RationalFunction(n, d)
Base.://(n::S, d::Poly{T}) where {T,S<:Number} = RationalFunction(n, d)
Base.://(n::Poly{T}, d::S) where {T,S<:Number} = RationalFunction(n, d)
# function Base.://(x::RationalFunction{T}, y::RationalFunction{T}) where T
#     xn, yn = divgcd(x.num, y.num)
#     xd, yd = divgcd(x.den, y.den)
#     (xn * yd) // (xd * yn)
# end
function Base.://(x::RationalFunction{T}, y::T) where T
    xn, yn = divgcd(x.num, y)
    xn // (x.den * yn)
end
function Base.://(x::T, y::RationalFunction{T}) where T
    xn,yn = divgcd(x, y.num)
    (xn * y.den) // yn
end

function Base.:*(x::RationalFunction, y::RationalFunction)
    x, y = promote(x, y)
    xn,yd = divgcd(x.num, y.den)
    xd,yn = divgcd(x.den, y.num)
    (xn * yn) // (xd * yd)
end

Base.:*(x::RationalFunction, y::Poly) = x * RationalFunction(y)
Base.:*(x::Poly{T}, y::RationalFunction{T}) where T = y * x
Base.:*(x::RationalFunction{T}, y::Number) where T = x * Poly(y, x.num.var)
Base.:*(x::Number, y::RationalFunction{T}) where T = y * x

# RationaFunction{T}(x::RationaFunction) where {T} = Rational{T}(convert(T,x.num), convert(T,x.den))
# RationaFunction{T}(x::Integer) where {T<:Integer} = Rational{T}(convert(T,x), convert(T,1))
# RationalFunction{T}(x::Poly{T}) = Rational{T}(convert(T,x), convert(T,1))

ispoly(x::RationalFunction) = x.den == 1

Base.numerator(x::RationalFunction) = x.num
Base.denominator(x::RationalFunction) = x.den

Base.promote_rule(::Type{RationalFunction{T}}, ::Type{RationalFunction{S}}) where {T,S} = RationalFunction{promote_type(S,T)}
Base.promote_rule(::Type{Poly{T}}, ::Type{T}) where {T} = Poly{T}

Base.convert(::Type{RationalFunction{T}}, x::RationalFunction{S}) where {T,S} = RationalFunction(convert(Poly{T}, numerator(x)), convert(Poly{T}, denominator(x)))
Base.convert(::Type{RationalFunction{T}}, x::Poly{T}) where {T} = RationalFunction(x)
Base.convert(::Type{RationalFunction{T}}, x::Number) where {T} = RationalFunction(T(x), T(1))

Base.convert(::Type{T}, r::RationalFunction{T}) where {T} = convert(T, numerator(r)) / convert(T, denominator(r))
Base.convert(::Type{T}, p::Poly{T}) where {T} = sum(c * Basic(p.var)^(i-1) for (i, c) in enumerate(coeffs(p)))

function monic(x::RationalFunction)
    n, d = numerator(x), denominator(x)
    lcn, lcd = lc(n), lc(d)
    lcn/lcd, n/lcn, d/lcd
end

shift(r::RationalFunction{T}, s::Union{Int64, T}) where {T} = RationalFunction(shift(numerator(r), s), shift(denominator(r), s))

(r::RationalFunction{T})(n::Union{Int, T}) where {T} = polyval(numerator(r), n) / polyval(denominator(r), n)

Base.one(::Type{RationalFunction{T}}) where {T} = one(T) // one(T)

# function convert(::Type{T}, r::RationalFunction{T}) where {T <: SymPy.Sym}
#     convert(T, numerator(r)) / convert(T, denominator(r))
# end