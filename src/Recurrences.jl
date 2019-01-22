module Recurrences

export test

using SymPy
using LinearAlgebra


include("zuercher.jl")
include("types.jl")

greet() = print("Hello World!")

function test()

    @vars n f g h k r v

    # cfhom = CFiniteRecurrence{Sym}(f, n, Sym[1, 2, 3])
    # cfinhom = CFiniteRecurrence{Sym}(f, n, Sym[3, 2, 1], Sym(4))

    # println(cfhom)
    # println(cfinhom)

    d1 = [Dict(f => Sym(2)), Dict(f => Sym(3))]
    d2 = [Dict(g => Sym(1)), Dict(f => Sym(-5))]
    d3 = [Dict(h => Sym(-1)), Dict(h => Sym(9))]
    d4 = [Dict(k => Sym(-1)), Dict(k => Sym(9)),Dict(k => Sym(-7))]
    
    lrs = LinearRecSystem(n)
    e1 = (coeffs = d1, inhom = Sym(2))
    e2 = (coeffs = d2, inhom = Sym(0))
    e3 = (coeffs = d3, inhom = Sym(-1))
    e4 = (coeffs = d4, inhom = Sym(32))
    push!(lrs, e1)
    push!(lrs, e2)
    push!(lrs, e3)
    push!(lrs, e4)

    display(lrs)

    d1 = [Dict(r => Sym(-1), v => Sym(1)), Dict(r => Sym(1))]
    d2 = [Dict(v => Sym(-1)), Dict(v => Sym(1))]
    
    lrs = LinearRecSystem(n)
    e1 = (coeffs = d1, inhom = Sym(0))
    e2 = (coeffs = d2, inhom = Sym(2))
    push!(lrs, e1)
    push!(lrs, e2)

    display(lrs)
    display(firstorder(lrs))
    display(homogenize!(lrs))
    display(decouple(lrs))
end

end # module
