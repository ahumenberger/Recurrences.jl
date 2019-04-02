import Base: convert, gcd
import SymPy.simplify

# ------------------------------------------------------------------------------

struct FallingFactorial{T}
    p::Poly{T}
    n::T
end

FallingFactorial(p::Poly{T}) where {T} = FallingFactorial(p, T(p.var))

shift(f::FallingFactorial{T}, s::Union{Int64, T}) where {T} = FallingFactorial(shift(f.p, s))
(f::FallingFactorial{T})(n::Union{Int, T}) where {T} = n == 0 ? 1 : prod(shift(f.p, -i)(n) for i in 0:n-1)

function convert(::Type{T}, f::FallingFactorial{T}) where {T <: SymPy.Sym}
    SymPy.FallingFactorial(convert(T, f.p), SymPy.Sym(f.p.var))
end

const RationalFactorial{T} = Pair{FallingFactorial{T},FallingFactorial{T}}
Base.numerator(f::RationalFactorial) = first(f)
Base.denominator(f::RationalFactorial) = last(f)

Base.one(::Type{RationalFactorial{T}}) where {T} = one(FallingFactorial{T}) => one(FallingFactorial{T})

Base.one(::Type{FallingFactorial{T}}) where {T} = FallingFactorial(one(Poly{T}), one(T))

convert(::Type{Basic}, f::FallingFactorial{Basic}) = Basic("gamma($(sprint(printpoly, f.p)) + 1) / gamma($(sprint(printpoly, f.p)) - $(f.n) + 1)")
convert(::Type{T}, f::RationalFactorial{T}) where {T} = convert(T, numerator(f)) / convert(T, denominator(f))

# ------------------------------------------------------------------------------

function convert(::Type{SymPy.Sym}, p::Poly)
    x = SymPy.Sym(string(p.var))
    sum(c*x^(i-1) for (i, c) in enumerate(coeffs(p)))
    # SymPy.Poly(sum(c*x^(i-1) for (i, c) in enumerate(coeffs(p))), x)
end

convert(::Type{Poly{SymPy.Sym}}, p::Poly{SymEngine.Basic}) = Poly(convert.(SymPy.Sym, coeffs(p)), p.var)
convert(::Type{Poly{SymEngine.Basic}}, p::Poly{SymPy.Sym}) = Poly(convert.(SymEngine.Basic, coeffs(p)), p.var)

convert(::Type{SymPy.Sym}, x::SymEngine.Basic) = begin 
    # @info "" x SymPy.Sym(string(x))
    SymPy.Sym(string(x))
end
convert(::Type{SymEngine.Basic}, x::SymPy.Sym) = begin
    # @info "asdf" x  convert(Expr, x)
    convert(Expr, x) |> SymEngine.Basic
end

convert(::Type{Poly}, p::SymPy.Sym) = Poly(SymPy.coeffs(p))

function resultant(p::Poly{T}, q::Poly{T}) where {T}
    res = SymPy.resultant(convert(SymPy.Sym, convert(Poly{Sym}, p)), convert(SymPy.Sym, convert(Poly{Sym}, q)))
    cf = coeffs(convert(T, res), convert(T, p.var))
    Poly(cf, p.var)
end

function gcd(p::Poly{SymPy.Sym}, q::Poly{SymPy.Sym})
    res = SymPy.gcd(convert(SymPy.Sym, p), convert(SymPy.Sym, q))
    if !SymPy.has(res, SymPy.Sym(p.var))
        return Poly([res], p.var)
    end
    Poly(SymPy.coeffs(SymPy.Poly(res, SymPy.Sym(p.var))), p.var)
end
gcd(p::Poly{SymEngine.Basic}, q::Poly{SymEngine.Basic}) = convert(Poly{Basic}, gcd(convert(Poly{Sym}, p), convert(Poly{Sym}, q)))

coeff(p::Poly{T}, x::T) where {T} = Poly(coeff.(coeffs(p), x), p.var)

Base.copysign(s::SymPy.Sym, i::Int64) = copysign(s, SymPy.Sym(i))
Base.copysign(s::SymPy.Sym, f::Float64) = copysign(s, SymPy.Sym(f))

function fallingfactorial(x, j)
    result = Poly([1], string(x))
    for i in 0:j - 1
        result *= Poly([-i, 1], string(x))
    end
    return result
end

SymPy.degree(p::Sym, n::Sym) = SymPy.degree(Poly(p, n))

function symset(v::String, j::Int64)
    return Sym[Sym("$v$i") for i in 1:j]
end

Base.denominator(s::SymPy.Sym) = denom(s)

function simplify(p::Poly)
    l = lcm2(denominator.(coeffs(p))...)
    Poly(coeffs(p) .* l, p.var)
end

function factors(expr::Sym)
    # c, list = factor_list(expr)
    # factor_list returns a Sym object instead of a tuple, does not seem to be right!
    # TOOO: add multiplicities
    c, list = factor_list(expr).x
    result = [x for (x,y) in list]
    # if c != 1
    #     push!(result, c)
    # end
    @debug "result" result
    return result
end

function shift(p::Poly{T}, s::Union{Int64, T}) where {T}
    c = T[polyval(polyder(p, i), s) / factorial(T(i)) for i in 0:degree(p)]
    Poly(c, p.var)
end

function factors(p::Poly{T}) where {T}
    sympoly = Sym(sprint(printpoly, p))
    facts = factors(sympoly)
    facts = convert(Vector{T}, facts)
    arg = convert(T, p.var)
    Poly.(coeffs.(facts, arg), p.var)
end

lc(p::Poly) = coeffs(p)[end]
coeff2(p::Poly, s::Int) = s <= degree(p) ? coeffs(p)[s+1] : 0

function mroots(poly::Poly{T}) where {T}
    roots = Polynomials.roots(poly)
    Dict([(r, count(x -> x==r, roots)) for r in Base.unique(roots)])
end

Base.promote_rule(::Type{T}, ::Type{SymEngine.Basic}) where{T<:SymPy.SymbolicObject} = SymEngine.Basic

Base.isequal(x::SymEngine.Basic, y::SymPy.Sym) = Base.isequal(promote(x, y)...)

mroots(p::Poly{SymPy.Sym}) = SymPy.polyroots(convert(SymPy.Sym, p))
mroots(p::Poly{SymEngine.Basic}) = convert(Dict{Basic,Int}, mroots(convert(Poly{SymPy.Sym}, p)))


# simplify(x::SymEngine.Basic) = SymEngine.expand(x)
simplify(x::SymEngine.Basic) = convert(SymEngine.Basic, SymPy.simplify(convert(SymPy.Sym, x)))

# Base.promote_rule(::Type{Basic}, ::Type{Int}) = Basic

# function mroots(p::Poly{SymPy.Sym})
    # T = SymPy.Sym
    # length(p) == 0 && return zeros(T, 0)

    # num_leading_zeros = 0
    # while p[num_leading_zeros] == zero(T)
    #     if num_leading_zeros == length(p)-1
    #         return zeros(T, 0)
    #     end
    #     num_leading_zeros += 1
    # end
    # num_trailing_zeros = 0
    # while p[end - num_trailing_zeros] == zero(T)
    #     num_trailing_zeros += 1
    # end
    # n = lastindex(p)-(num_leading_zeros + num_trailing_zeros)
    # n < 1 && return zeros(T, length(p) - num_trailing_zeros - 1)

    # companion = diagm(-1 => ones(T, n-1))
    # an = p[end-num_trailing_zeros]
    # companion[1,:] = -p[(end-num_trailing_zeros-1):-1:num_leading_zeros] / an

    # companion[:eigenvals]()
# end

function linsolve(A::Matrix{Basic}, b::Vector{Basic})
    nrows = size(A, 1)
    ncols = size(A, 2)
    for i in nrows+1:ncols
        v = zeros(Basic, 1, ncols)
        v[1, i] = 1
        A = vcat(A, v)
        b = [b; 0]
    end
    cA = convert(SymEngine.CDenseMatrix, A)
    # @info "" A cA b
    res = convert(Matrix{Basic}, cA \ b)
    vec(res)
end