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
include("parse.jl")
include("sympy.jl")
include("utils.jl")
include("cfs.jl")

const AppliedUndef = PyCall.PyNULL()

function __init__()
    copy!(AppliedUndef, PyCall.pyimport_conda("sympy.core.function", "sympy").AppliedUndef)
end

end # module
