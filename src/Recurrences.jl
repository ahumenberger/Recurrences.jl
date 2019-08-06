module Recurrences

export @rec, @lrs, @lrs_parallel, @solve_parallel
export lrs, lrs_sequential, lrs_parallel
export Recurrence, LinearRecurrence, CFiniteRecurrence, HyperRecurrence
export ClosedForm, CFiniteClosedForm, HyperClosedForm
export closedform, expression, exponentials
export LinearRecSystem, decouple, homogenize!, solve
export initvar

using SymPy
using SymEngine
using PyCall
using LinearAlgebra
using Polynomials
using MacroTools
using DynamicPolynomials
using MultivariatePolynomials

import Polynomials: Poly, printpoly, degree, coeffs, polyval, polyder
import Base: convert, denominator

const RExpr = Union{Expr,Symbol,Number}
const APL = AbstractPolynomialLike
const Var = AbstractVariable
const RPoly = Polynomial{true,Rational{Int}}

include("rationalfunction.jl")
include("polyhelpers.jl")
include("zuercher.jl")
include("petkovsek.jl")
include("rectypes.jl")
include("cftypes.jl")
include("recsystem.jl")
include("sympy.jl")
include("utils.jl")
include("poly.jl")

const AppliedUndef = PyCall.PyNULL()

function __init__()
    copy!(AppliedUndef, PyCall.pyimport_conda("sympy.core.function", "sympy").AppliedUndef)
end

replace_post(ex, s, s′) = MacroTools.postwalk(x -> x == s ? s′ : x, ex)
gensym_unhashed(s::Symbol) = Symbol(replace(string(gensym(s)), "#"=>""))

function split_assign(xs::Vector{Expr})
    ls = Symbol[]
    rs = RExpr[]
    for x in xs
        @capture(x, l_ = r_)
        push!(ls, unblock(l))
        push!(rs, unblock(r))
    end
    ls, rs
end

function pushexpr!(lrs::LinearRecSystem, ls::Vector{Expr}, rs::AbstractVector{RExpr})
    for (l, r) in zip(ls, rs)
        expr = :($l - $r)
        entry, _ = LinearRecEntry(Basic, expr)
        push!(lrs, entry)
    end
    lrs
end

function pushexpr!(lrs::LinearRecSystem{S,T}, ex::Expr...) where {S,T}
    for x in ex
        if @capture(x, l_ = r_)
            x = :($(unblock(l)) - $(unblock(r)))
        end
        entry, _ = LinearRecEntry(Var, APL, x)
        @info "" entry.coeffs entry.inhom
        @info "" entry isa LinearEntry{S,T}
        # @info "" LinearEntry{S,T}(entry) lrs
        push!(lrs, entry)
    end
    lrs
end

function _parallel(lhss::Vector{Symbol}, rhss::Vector{RExpr})
    del = Int[]
    for (i,v) in enumerate(lhss)
        j = findnext(x->x==v, lhss, i+1)
        l = j === nothing ? lastindex(lhss) : j
        for k in i+1:l
            rhss[k] = replace_post(rhss[k], v, rhss[i])
        end
        if j != nothing
            push!(del, i)
        end
    end
    
    if !isempty(del)
        deleteat!(lhss, Tuple(del))
        deleteat!(rhss, Tuple(del))
    end
    lhss, rhss
end

function lrs_sequential(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lhss, rhss = _parallel(split_assign(exprs)...)
    @info "Splitted and parallel assignments" rhss lhss
    _lrs_parallel(lhss, rhss, lc)
    # lrs = LinearRecSystem(Basic(lc), map(Basic, Base.unique(lhss)))
    # for (i,v) in enumerate(lhss)
    #     rhss = RExpr[replace_post(x, v, (i < j ? :($v($lc+1)) : :($v($lc)))) for (j,x) in enumerate(rhss)]
    # end
    # lhss = [Expr(:call, v, Expr(:call, :+, lc, 1)) for v in lhss]
    # pushexpr!(lrs, lhss, rhss)
end

function lrs_parallel(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lhss, rhss = split_assign(exprs)
    _lrs_parallel(lhss, rhss, lc)
end

function _lrs_parallel(lhss::AbstractVector{Symbol}, rhss::AbstractVector{RExpr}, lc::Symbol)
    linear = findall(x->_islinear(x, lhss), rhss)
    nonlinear = setdiff(eachindex(lhss), linear)
    rest = (view(lhss, nonlinear), view(rhss, nonlinear))
    _lhss, _rhss = view(lhss, linear), view(rhss, linear)
    lrs = LinearRecSystem(Basic(lc), map(Basic, _lhss))
    for (i, rhs) in enumerate(_rhss)
        _rhss[i] = MacroTools.postwalk(x -> x isa Symbol && x in _lhss ? :($x($lc)) : x, rhs)
    end
    _lhss = [Expr(:call, v, Expr(:call, :+, lc, 1)) for v in _lhss]
    pushexpr!(lrs, _lhss, _rhss), rest
end

function solve_sequential(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))

end

function solve_parallel(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lhss, rhss = split_assign(exprs)
    _solve_parallel(lhss, rhss, lc)
end

function _solve_parallel(lhss::AbstractVector{Symbol}, rhss::AbstractVector{RExpr}, lc::Symbol)
    lrs, remaining = _lrs_parallel(lhss, rhss, lc)
    @info "" lrs
    cfs = solve(lrs)
    cfmap = Dict{Symbol, RExpr}(Symbol(string(x.func)) => convert(Expr, expression(x)) for x in cfs)
    _lhss, _rhss = remaining
    isempty(_lhss) && return (cfs, remaining)

    for (i, rhs) in enumerate(_rhss)
        _rhss[i] = MacroTools.postwalk(x -> x isa Symbol && x in keys(cfmap) ? cfmap[x] : x, rhs)
    end
    _cfs, _ = _solve_parallel(_lhss, _rhss, lc)
    if isempty(_cfs)
        @warn("Cannot solve $(_lhss)")
    end
    return cfs, remaining
end

function lrs(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    @info "" exprs lc
    lrs = LinearRecSystem{Var,APL}(mkvar(lc))
    pushexpr!(lrs, exprs...)
end

end # module
