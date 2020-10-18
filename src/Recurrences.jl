module Recurrences

export @lrs
export Seq, ClosedForm

using MacroTools: postwalk, @capture, unblock
using Combinatorics
import AbstractAlgebra: Ring, RingElement, Field, FieldElement, Generic.Frac, RingElem
import AbstractAlgebra: parent_type, elem_type, base_ring, parent, needs_parentheses, displayed_with_minus_in_front, canonical_unit
import AbstractAlgebra: zero!, mul!, add!, addeq!
using Nemo
using SymEngine
using LinearAlgebra
import Base: convert

const RExpr = Union{Expr,Symbol,Number}

include("helper.jl")
include("utils.jl")
include("hyperterm.jl")
include("sequence.jl")
include("zuercher.jl")
include("rectypes.jl")
include("recsystem.jl")
include("petkovsek.jl")
include("macros.jl")

end # module