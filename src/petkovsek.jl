# module Petkovsek

# using SymPy

# export algpoly

function fallingfactorial(x, j)
    result = Polynomials.Poly([1])
    for i in 0:j - 1
        result *= Polynomials.Poly([-i, 1])
    end
    return result
end

SymPy.degree(p::Sym, n::Sym) = SymPy.degree(Poly(p, n))

function symset(v::String, j::Int64)
    return Sym[Sym("$v$i") for i in 1:j]
end



function factors(expr::Sym)
    # c, list = factor_list(expr)
    # factor_list returns a Sym object instead of a tuple, does not seem to be right!
    c, list = factor_list(expr).x
    result = [x^y for (x,y) in list]
    # if c != 1
        # push!(result, c)
    # end
    return result
end

function shift(p::Polynomials.Poly{T}, s::Int) where{T}
    c = T[polyval(polyder(p, i), s) / factorial(i) for i in 0:Polynomials.degree(p)]
    Polynomials.Poly(c, p.var)
end

function factors(p::Polynomials.Poly)
    sympoly = Sym(sprint(printpoly, p))
    facts = factors(sympoly)
    Polynomials.Poly.(SymPy.coeffs.(facts))
end

lc(p::Polynomials.Poly) = Polynomials.coeffs(p)[end]
coeff2(p::Polynomials.Poly, s::Int) = s <= Polynomials.degree(p) ? Polynomials.coeffs(p)[s+1] : 0


function algpoly(polylist::Vector{Polynomials.Poly{T}}, f, n) where {T}
    @debug "algpoly" polylist
    # Construct polynomials qj that arise from treating the shift operator
	# as the difference operator plus the identity operator.
    qlist = []
    for j in 0:length(polylist)-1
        qj = 0 * n
        for i in j:length(polylist)-1
            qj += binomial(i, j) * polylist[i+1]
        end
        push!(qlist, qj)
    end
    @debug "haha" [(Polynomials.degree(poly)) for (j, poly) in enumerate(qlist)]
    # Find all candidates for the degree bound on output polynomial.
    b = maximum([(Polynomials.degree(poly) - (j - 1)) for (j, poly) in enumerate(qlist)])
    @debug "" qlist b
    first = Polynomials.degree(f) - b
    second = -1 * b - 1
    lcoeffs = lc.(qlist)
    @debug "" lcoeffs
    alpha = 0
    for j in 0:length(qlist)-1
        aux = Polynomials.degree(qlist[j+1]) - j
        @debug "what" aux b aux == b
        if (Polynomials.degree(qlist[j+1]) - j) == b
            alpha += lcoeffs[j+1] * fallingfactorial(n, j)
        end
    end
    @debug "" alpha
    deg = mroots(alpha)
    third = maximum(keys(deg))
    d = convert(Int64, max(first, second, third, 0))
    @debug "Step 2" d
    # Use method of undetermined coefficients to find the output polynomial.
    varlist = symset("a", d+1)
    pol = 0 * n
    for i in 0:d
        pol += n^i * varlist[i+1]
    end
    poly = Polynomials.Poly(varlist)
    # solution = -1 * f
    # for (i, p) in enumerate(polylist)
    #     @debug "" shift(poly, i - 1) * p
    #     # solution += subs(pol, (n, n + (i-1))) * poly
    #     # @debug "" solution
    # end

    sh = [shift(poly, i - 1) * p for (i, p) in enumerate(polylist)]
    maxdeg = maximum(Polynomials.degree(p) for p in sh)
    A = [coeff2(p, d) for d in 0:maxdeg, p in sh]
    b = [coeff2(p, 0) for p in sh]
    @debug "" sh A b #A\b
        # solution += subs(pol, (n, n + (i-1))) * poly
        # @debug "" solution

    # coef = Polynomials.coeffs(solution)
    # filter!(e->e!=n*0, coef)
    # if isempty(coef)
    #     return Dict([v => v for v in varlist])
    # end
    # return solve(coef, varlist)
    # A \ b
end

# function Base.convert(::Type{SymPy.Sym}, p::Polynomials.Poly{SymPy.Sym})

# end

Base.copysign(s::SymPy.Sym, i::Int64) = copysign(s, SymPy.Sym(i))
Base.copysign(s::SymPy.Sym, f::Float64) = copysign(s, SymPy.Sym(f))


function alghyper(polylist::Vector{Polynomials.Poly{T}}, n::T) where {T}
    p = polylist[1]
    alist = factors(p)
    lcoeff = lc(p)
    for (i,p) in enumerate(alist)
        alist[i] /= lcoeff
    end
    if !(Polynomials.Poly([T(1)]) in alist)
        push!(alist, Polynomials.Poly([T(1)]))
    end

    d = length(polylist)
    p = shift(polylist[end], -d+2)
    lcoeff = lc(p)
    blist = factors(p)
    for (i,p) in enumerate(blist)
        blist[i] /= lcoeff
    end
    if !(Polynomials.Poly([T(1)]) in blist)
        push!(blist, Polynomials.Poly([T(1)]))
    end

    solutions = []
    for aelem in alist
        for belem in blist
            plist = []
            for i in 0:d-1 
                pi = polylist[i+1]
                for j in 0:i-1
                    pi *= shift(aelem, j)
                end
                for j in i:d-1
                    pi *= shift(belem, j)
                end
                push!(plist, pi)
            end

            m = maximum(Polynomials.degree.(plist))
            alpha = coeff2.(plist, m)
            @syms z
            zpol = 0*z
            for i in 0:length(alpha) - 1
                zpol += alpha[i+1]*z^i
            end

            vals = [key for (key,val) in polyroots(zpol) if key != 0]
            @debug "roots" vals
            for x in vals
                polylist2 = [x^(i-1)*p for (i,p) in enumerate(plist)]
                
                polysols = algpoly(polylist2, Polynomials.Poly([0]), n)
                if isempty(polysols)
                    continue
                end
                # polysols = collect(values(polysols))
                filter!(e->e!=n*0, polysols)
                @debug "" polysols
                if length(polysols) > 0
                    c = 0*n
                    for (i,p) in enumerate(polysols)
                        c += n^(i-1) * p
                    end
                    s = x * aelem/belem * subs(c, (n, n+1))/c
                    push!(solutions, simplify(s))
                end
            end
        end
    end
    return solutions
end

function tohg(sol, n)
    facts = factors(sol)
    result = []
    for f in facts
        if has(f, n)
            f = f |> subs(n, n-1)
            c = coeff(f, n)
            f = f / c
            push!(result, factorial(f))
            if c != 1
                push!(result, c^n)
            end
        else
            push!(result, f^n)
        end
    end
    return prod(result)
end

function summands(expr::Sym)
    expr = expand(expr)
    @debug "" expr typeof(expr)
    if SymPy.funcname(expr) == "Add"
        return args(expr)
    end
    [expr]
end

lcm2(n::Sym, rest::Sym...) = SymPy.lcm(n, lcm2(rest...))
lcm2() = 1

function petkovsek(coeffs::Vector{T}, arg::T) where {T}
    @debug "" coeffs
    ls = summands.(coeffs) |> Iterators.flatten
    @debug "" ls |> collect
    ds = denom.(ls)
    val = lcm2(ds...)
    coeffs *= val
    coeffs = simplify.(coeffs)
    coeffs = SymPy.Poly.(coeffs, arg)
    @debug "" coeffs
    hyper = alghyper(Polynomials.Poly.(SymPy.coeffs.(coeffs)), arg)
    @debug "Petkovsek - alghyper" hyper
    tohg.(hyper, arg)
end

# function hyper()

# end # module

# using SymPy
# using Petkovsek

# @syms n y

# alghyper([2*n*(n+1), -(n^2 +3*n-2), n-1], n)

# # algpoly([3, -n, n-1], 0*n, n)
# # algpoly([n-1, -n, 3], 0*n, n)
# println("alghyper: ", [tohg(sol, n) for sol in alghyper([2*n*(n+1), -(n^2 +3*n-2), n-1], n)])
# println("algpoly: ", algpoly([n*(n + 1), -1*n^2 - 3*n + 2, 2*n - 2], 0*n,n))