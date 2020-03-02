module Recurrences

export @lrs

using MacroTools: postwalk, @capture, unblock
using Combinatorics
import AbstractAlgebra: Ring, RingElement, Field, FieldElement
import AbstractAlgebra: parent_type, elem_type, base_ring, parent, needs_parentheses, displayed_with_minus_in_front, canonical_unit
using Nemo
using SymEngine
using LinearAlgebra
import Base: convert

const RExpr = Union{Expr,Symbol,Number}

include("new/helper.jl")
include("new/utils.jl")
include("new/zuercher.jl")
include("new/rectypes.jl")
include("new/recsystem.jl")
include("new/hyperterm.jl")
include("new/sequence.jl")
include("new/petkovsek.jl")
include("new/macros.jl")

end # module