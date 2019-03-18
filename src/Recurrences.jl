module Recurrences

export lrs_sequential

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

function lrs_sequential(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lhss = Symbol[]
    rhss = Expr[]
    for assign in exprs
        @capture(assign, lhs_ = rhs_)
        push!(lhss, lhs)
        push!(rhss, unblock(rhs))
    end

    sys = LinearRecSystem(Basic(lc), map(Basic, lhss))

    for (i,v) in enumerate(lhss)
        rhss = [replace_post(x, v, (i < j ? :($v($lc+1)) : :($v($lc)))) for (j,x) in enumerate(rhss)]
    end
    lhss = [Expr(:call, v, Expr(:call, :+, lc, 1)) for v in lhss]

    entries = Recurrences.LinearEntry{Basic}[]    
    for (lhs, rhs) in zip(lhss, rhss)
        expr = Expr(:call, :(-), lhs, rhs)
        entry, _ = Recurrences.LinearRecEntry(Basic, expr)
        push!(entries, entry)
    end
    push!(sys, entries...)
    @debug "Rec system" sys
    sys
end

end # module
