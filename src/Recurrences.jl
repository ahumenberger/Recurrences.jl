module Recurrences

export @rec, @lrs
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

include("rationalfunction.jl")
include("polyhelpers.jl")
include("zuercher.jl")
include("petkovsek.jl")
include("rectypes.jl")
include("cftypes.jl")
include("recsystem.jl")
include("sympy.jl")
include("utils.jl")

const AppliedUndef = PyCall.PyNULL()

function __init__()
    copy!(AppliedUndef, PyCall.pyimport_conda("sympy.core.function", "sympy").AppliedUndef)
end

replace_post(ex, s, s′) = MacroTools.postwalk(x -> x == s ? s′ : x, ex)
gensym_unhashed(s::Symbol) = Symbol(replace(string(gensym(s)), "#"=>""))

const RExpr = Union{Expr,Symbol,Number}

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

function pushexpr!(lrs::LinearRecSystem, ls::Vector{Expr}, rs::Vector{RExpr})
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
        deleteat!(lhss, del...)
        deleteat!(rhss, del...)
    end
    lhss, rhss
end

# _islinear(expr, vars) = SymEngine.degree() degree()

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

function _lrs_parallel(lhss::Vector{Symbol}, rhss::Vector{RExpr}, lc::Symbol)
    lrs = LinearRecSystem(Basic(lc), map(Basic, lhss))
    for (i, rhs) in enumerate(rhss)
        rhss[i] = MacroTools.postwalk(x -> x isa Symbol && x in lhss ? :($x($lc)) : x, rhs)
    end
    lhss = [Expr(:call, v, Expr(:call, :+, lc, 1)) for v in lhss]
    pushexpr!(lrs, lhss, rhss)
end

function lrs(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lrs = LinearRecSystem(Basic(lc))
    pushexpr!(lrs, exprs...)
end

end # module
