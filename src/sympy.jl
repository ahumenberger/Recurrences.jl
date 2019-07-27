free_symbols(x::Union{SymPy.Sym, Vector{SymPy.Sym}}) = SymPy.free_symbols(x)
free_symbols(x::Union{SymEngine.Basic, Vector{SymEngine.Basic}}) = SymEngine.free_symbols(x)

coeff(x::SymPy.Sym, b::SymPy.Sym) = x.coeff(b)

subs(x::SymPy.Sym, args...) = SymPy.subs(x, args...)
subs(x::SymEngine.Basic, args...) = SymEngine.subs(x, args...)

subs(ex::AbstractArray, args...; kwargs...) = map(u -> subs(u, args...; kwargs...), ex)
# subs(ex::AbstractArray{SymEngine.Basic}, args...; kwargs...) = map(u -> subs(u, args...; kwargs...), ex)

Base.inv(m::Matrix{SymEngine.Basic}) = convert(Matrix{SymEngine.Basic}, inv(convert(Matrix{Sym}, m)))

expand(x::Sym) = SymPy.expand(x)
expand(x::Basic) = SymEngine.expand(x)

solve(x::Vector{Sym}, y::Vector{Sym}) = SymPy.solve(x, y)
solve(x::Vector{Basic}, y::Vector{Basic}) = convert(Dict{Basic,Basic}, solve(convert.(Sym, x), convert.(Sym, y)))

coeff(x::Basic, y::Basic) = SymEngine.coeff(x, y)
coeff(x::Basic, y::Basic) = convert(Basic, coeff(convert(Sym, x), convert(Sym, y)))

function coeffs(p::Sym, x::Sym)
    if iszero(p)
        return []
    end
    SymPy.coeffs(p, x)
end
function coeffs(p::Basic, x::Basic)
    convert.(Basic, coeffs(convert(Sym, p), convert(Sym, x)))
end

Base.:\(A::Matrix{Basic}, b::Vector{Basic}) = convert(Vector{Basic}, convert(Matrix{Sym}, A) \ convert(Vector{Sym}, b))