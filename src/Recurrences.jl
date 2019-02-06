module Recurrences

export test, test2, test3

using SymPy
using PyCall
using LinearAlgebra
using Polynomials
using MacroTools

import Polynomials: Poly, printpoly, degree, coeffs, polyval, polyder

include("rationalfunction.jl")
include("zuercher.jl")
include("petkovsek.jl")
include("rectypes.jl")
include("cftypes.jl")
include("recsystem.jl")
include("sympy.jl")
include("utils.jl")

greet() = print("Hello World!")

const AppliedUndef = PyCall.PyNULL()

function __init__()
    copy!(AppliedUndef, PyCall.pyimport_conda("sympy.core.function", "sympy")[:AppliedUndef])
end

function test()

    @vars n f g h k r v

    # cfhom = CFiniteRecurrence{Sym}(f, n, Sym[1, 2, 3])
    # cfinhom = CFiniteRecurrence{Sym}(f, n, Sym[3, 2, 1], Sym(4))

    # println(cfhom)
    # println(cfinhom)

    d1 = [Dict(f => Sym(2)), Dict(f => Sym(3))]
    # d2 = [Dict(g => Sym(1)), Dict(f => Sym(-5))]
    # d3 = [Dict(h => Sym(-1)), Dict(h => Sym(9))]
    # d4 = [Dict(k => Sym(-1)), Dict(k => Sym(9)),Dict(k => Sym(-7))]
    
    # lrs = LinearRecSystem(n)
    # e1 = (coeffs = d1, inhom = Sym(2))
    # e2 = (coeffs = d2, inhom = Sym(0))
    # e3 = (coeffs = d3, inhom = Sym(-1))
    # e4 = (coeffs = d4, inhom = Sym(32))
    # push!(lrs, e1)
    # push!(lrs, e2)
    # push!(lrs, e3)
    # push!(lrs, e4)

    # display(lrs)

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
    lrs = decouple(lrs)
    display(lrs)


    # firstorder(lrs)
    # homogenize!(lrs)
    # decouple(lrs)
end

function test2()
    @vars a b c d e
    @vars x y
    @vars n

    # q = a*c*e/(d*(b*d + c*e - e*(a + c) + e)) + c/d + e*(-a - c)/(d*(b*d + c*e - e*(a + c) + e)) - 1/d*x + e/(d*(b*d + c*e - e*(a + c) + e))*x^2
    # q = simplify(q)
    # p = a*c + -a*c - a - c*x + a + c + 1*x^2 + -1*x^3

    # @info "" p q
    # q = -a*e*x + a*e - b*c*d + b*d*x - e*x^2 + e*x
    # @info "" q
    # @info "" divrem(Polynomials.Poly(SymPy.coeffs(p)), Polynomials.Poly(SymPy.coeffs(q)))

    # @info "" simplify(c + (1 - (-a*e + b*d + e)/e)*(a*e - b*c*d)/e + -c + (1 - (-a*e + b*d + e)/e)*(-a*e + b*d + e)/e - (a*e - b*c*d)/e*x)

    d1 = [Dict(x => Sym(-a)), Dict(x => Sym(1))]
    # e1 = (coeffs = d1, inhom = b)
    e1 = (coeffs = d1, inhom = Sym(0))
    d2 = [Dict(y => Sym(-c)), Dict(y => Sym(1), x => Sym(-d))]
    # e2 = (coeffs = d2, inhom = e)
    e2 = (coeffs = d2, inhom = Sym(0))
    @info "Linear entries" e1 e2

    lrs = LinearRecSystem(n)
    push!(lrs, e2)
    push!(lrs, e1)
    @info "Linear rec system" lrs

    lrs = decouple(lrs)
    @info "Decoupled" lrs
end

function test3()
    @vars a b c d e
    @vars x y
    @vars n



    d1 = [Dict(x => Sym(-2)), Dict(x => Sym(1))]
    # e1 = (coeffs = d1, inhom = b)
    e1 = (coeffs = d1, inhom = Sym(0))
    d2 = [Dict(y => Sym(-3), x => Sym(-4)), Dict(y => Sym(1))]
    # e2 = (coeffs = d2, inhom = e)
    e2 = (coeffs = d2, inhom = Sym(0))
    @info "Linear entries" e1 e2

    lrs = LinearRecSystem(n)
    push!(lrs, e2)
    push!(lrs, e1)
    @info "Linear rec system" lrs

    lrs = decouple(lrs)
    @info "Decoupled" lrs
end

end # module
