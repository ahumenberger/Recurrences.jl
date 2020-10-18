export mkvar, mkpoly

mkpoly(x::Expr) = eval(symbol_walk(y->mkvar(y), x))
mkpoly(x::String) = mkpoly(Meta.parse(x))
mkpoly(x::Symbol) = mkvar(x)
mkpoly(x::Number) = iszero(x) ? zero(Polynomial{true,Int}) : one(Polynomial{true,Int}) * x

_varmap = Dict{Symbol, AbstractVariable}()

mkvar(x::String) = mkvar(Symbol(x))
mkvar(x::T) where {T <: AbstractVariable} = x
function mkvar(x::Symbol)
    if haskey(_varmap, x)
        return _varmap[x]
    end
    v = PolyVar{true}(string(x))
    push!(_varmap, x=>v)
    v
end

function _islinear(x::RExpr, vars::AbstractVector{Symbol})
    vs = map(mkvar, vars)
    p = mkpoly(x)
    all(MultivariatePolynomials.degree(m, v) < 2 for v in vs for m in MultivariatePolynomials.terms(p))
end

# Base.zeros(::Type{APL}, dims::Base.DimOrInd...) = zeros(APL, dims...)
# Base.zero(::Type{APL}) = zero(APL)

# Base.promote_op(transpose, ::Type{APL}...) = APL
# Base.promote_op(transpose, ::Type{RPoly}...) = RPoly

# APL(x::Bool) = x ? one(APL) : zero(APL)

# LinearAlgebra.transpose(p::RationalPoly) = transpose(numerator(p)) / transpose(denominator(p))
# LinearAlgebra.adjoint(p::RationalPoly) = adjoint(numerator(p)) / adjoint(denominator(p))