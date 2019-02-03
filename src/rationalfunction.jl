export RationalFunction
# import Base: //

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
RationalFunction(n::Poly{T}) where {T} = RationalFunction(n, one(T))
RationalFunction(x::RationalFunction) = x

function divgcd(x::Poly{T}, y::Poly{T}) where T
    g = gcd(x, y)
    div(x, g), div(y, g)
end

Base.://(n::Poly{T},  d::Poly{T}) where T = RationalFunction(n, d)
Base.://(n::Poly{T},  d::S) where {T,S<:Number} = RationalFunction(n, d)
function Base.://(x::RationalFunction{T}, y::RationalFunction{T}) where T
    xn, yn = divgcd(x.num, y.num)
    xd, yd = divgcd(x.den, y.den)
    (xn * yd) // (xd * yn)
end
function Base.://(x::RationalFunction{T}, y::T) where T
    xn, yn = divgcd(x.num, y)
    xn // (x.den * yn)
end
function Base.://(x::T,  y::RationalFunction{T}) where T
    xn,yn = divgcd(x, y.num)
    (xn * y.den) // yn
end

function Base.:*(x::RationalFunction, y::RationalFunction)
    xn,yd = divgcd(x.num, y.den)
    xd,yn = divgcd(x.den, y.num)
    (xn * yn) // (xd * yd)
end

Base.:*(x::RationalFunction{T}, y::Poly{T}) where T = x * (y // 1)
Base.:*(x::Poly{T}, y::RationalFunction{T}) where T = y * x
Base.:*(x::RationalFunction{T}, y::T) where T = x * Poly(y, x.num.var)
Base.:*(x::T, y::RationalFunction{T}) where T = y * x

# RationaFunction{T}(x::RationaFunction) where {T} = Rational{T}(convert(T,x.num), convert(T,x.den))
# RationaFunction{T}(x::Integer) where {T<:Integer} = Rational{T}(convert(T,x), convert(T,1))
# RationalFunction{T}(x::Poly{T}) = Rational{T}(convert(T,x), convert(T,1))

ispoly(x::RationalFunction) = x.den == 1

numerator(x::RationalFunction) = x.num
denominator(x::RationalFunction) = x.den

Base.promote_rule(::Type{RationalFunction{T}}, ::Type{RationalFunction{S}}) where {T,S} = RationalFunction{promote_type(S,T)}
Base.promote_rule(::Type{Poly{T}}, ::Type{T}) where {T} = Poly{T}

Base.convert(::Type{RationalFunction{T}}, x::Poly{T}) where {T} = RationalFunction(x)

function monic(x::RationalFunction)
    n, d = numerator(x), denominator(x)
    lcn, lcd = lc(n), lc(d)
    lcn/lcd, n/lcn, d/lcd
end