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

import Polynomials: Poly, printpoly, degree, coeffs, polyval, polyder
import Base: convert, denominator

const RExpr = Union{Expr,Symbol,Number}

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
replace_post(ex, dict) = MacroTools.postwalk(x -> x in keys(dict) ? dict[x] : x, ex)
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

function pushexpr!(lrs::LinearRecSystem, ex::Expr...)
    for x in ex
        if @capture(x, l_ = r_)
            x = :($l - $r)
        end
        entry, _ = LinearRecEntry(Basic, x)
        push!(lrs, entry)
    end
    lrs
end

function _parallel(lhss::Vector{Symbol}, rhss::Vector{RExpr})
    del = Int[]
    dict = collect(zip(lhss, rhss))
    for (i,v) in enumerate(lhss)
        rhss[i] = replace_post(rhss[i], Dict(dict[1:i-1]))
    end
    _lhss, _rhss = Symbol[], RExpr[]
    for (l, r) in reverse(collect(zip(lhss, rhss)))
        if l ∉ _lhss
            pushfirst!(_lhss, l)
            pushfirst!(_rhss, r)
        end
    end
    _lhss, _rhss
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
    lrs = LinearRecSystem(Basic(lc))
    pushexpr!(lrs, exprs...)
end

end # module
