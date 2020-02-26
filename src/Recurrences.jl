module Recurrences

export @lrs

using MacroTools: postwalk, @capture, unblock
# using AbstractAlgebra
using Nemo
using SymEngine
using LinearAlgebra
import Base: convert

const RExpr = Union{Expr,Symbol,Number}

include("new/helper.jl")
include("new/utils.jl")
include("new/zuercher.jl")
include("new/petkovsek.jl")
include("new/rectypes.jl")
include("new/cftypes.jl")
include("new/recsystem.jl")
include("new/macros.jl")

end # module