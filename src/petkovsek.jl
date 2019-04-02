function algpoly(polylist::Vector{Poly{T}}, f::Poly{T}, n) where {T}
    @debug "algpoly" polylist
    # Construct polynomials qj that arise from treating the shift operator
	# as the difference operator plus the identity operator.
    qlist = []
    for j in 0:length(polylist)-1
        qj = 0
        for i in j:length(polylist)-1

            qj += (binomial(i, j) * polylist[i+1])

        end
        push!(qlist, qj)
    end

    # Find all candidates for the degree bound on output polynomial.
    b = maximum([(degree(poly) - (j - 1)) for (j, poly) in enumerate(qlist)])
    first = degree(f) - b
    second = -1 * b - 1
    lcoeffs = lc.(qlist)
    alpha = 0
    for j in 0:length(qlist)-1
        aux = degree(qlist[j+1]) - j
        if (degree(qlist[j+1]) - j) == b
            alpha += lcoeffs[j+1] * fallingfactorial(n, j)
        end
    end
    @debug "" alpha typeof(alpha)
    deg = mroots(alpha)
    deg = isa(deg, Dict) ? keys(deg) : deg
    third = convert(Int64, maximum(deg))
    @debug "" qlist b deg first second third
    d = max(first, second, third, 0)

    # Use method of undetermined coefficients to find the output polynomial.
    varlist = variables(T, n=d+1)
    # pol = 0 * n
    # for i in 0:d
    #     pol += n^i * varlist[i+1]
    # end
    genpoly = Poly(varlist, string(n))
    @debug "Generic polynomial" genpoly
    # solution = -1 * f
    # for (i, p) in enumerate(polylist)
    #     @debug "" shift(poly, i - 1) * p
    #     # solution += subs(pol, (n, n + (i-1))) * poly
    #     # @debug "" solution
    # end

    poly = sum(shift(genpoly, i - 1) * p for (i, p) in enumerate(polylist))
    poly -= f
    @debug "Poly" poly
    if iszero(poly)
        return Poly(one(T), string(n))
    end
    coefficients = coeffs(poly)
    filter!(e -> e != 0, coefficients)
    sol = solve(coefficients, varlist)
    missing = setdiff(varlist, keys(sol))
    p = Poly([subs(c, sol...) for c in coeffs(genpoly)], string(n))
    solutions = Poly{T}[]
    for v in missing
        c = coeff(p, v)
        if !iszero(c)
            push!(solutions, c)
        end
    end
    @debug "Solution for coefficients" sol p solutions
    return solutions
end

function alghyper(polylist::Vector{Poly{T}}, n::T) where {T}
    @debug "-> alghyper" polylist n
    p = polylist[1]
    alist = factors(p)
    lcoeff = lc(p)

    # for (i,p) in enumerate(alist)
    #     alist[i] /= lcoeff
    # end
    if !(1 in alist)
        push!(alist, Poly([T(1)], string(n)))
    end

    d = length(polylist)
    p = shift(polylist[end], -d+2)
    lcoeff = lc(p)
    blist = factors(p)
    # for (i,p) in enumerate(blist)
    #     blist[i] /= lcoeff
    # end
    if !(1 in blist)
        push!(blist, Poly([T(1)], string(n)))
    end

    alist, blist = reverse(alist), reverse(blist)

    @debug "" alist blist

    solutions = []
    for aelem in alist
        for belem in blist
            @debug "Monic factors" aelem belem
            plist = Poly{T}[]
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

            m = maximum(degree.(plist))
            alpha = convert(Vector{T}, coeff2.(plist, m))
            # @syms z
            # zpol = 0*z
            # for i in 0:length(alpha) - 1
            #     zpol += alpha[i+1]*z^i
            # end
            zpol = Poly(alpha)

            @debug "" zpol typeof(zpol) typeof(alpha) typeof(plist)
            vals = [key for (key,val) in mroots(zpol) if key != 0]
            @debug "Roots" zpol vals
            for x in vals
                polylist2 = [x^(i-1)*p for (i,p) in enumerate(plist)]
                @debug "Input to algpoly" polylist2
                
                polysols = algpoly(polylist2, Poly(zero(T), string(n)), n)
                @debug "Solutions of algpoly" polysols

                for c in polysols
                    @debug "" c
                    s = x * (aelem // belem) * (shift(c, 1) // c)
                    @info "Solution" s
                    push!(solutions, s)
                end
            end
        end
    end
    return solutions
end

function commonfactors(p::Poly{T}, q::Poly{T}) where {T}
    if p == 1 || q == 1
        return Poly(one(T), p.var) // 1, Pair(FallingFactorial(p), FallingFactorial(q))
    end
    k = variables(T)
    res = resultant(shift(p, k), q)
    roots = keys(mroots(res)) |> collect
    filter!(x -> isinteger(x), roots)
    @debug "Integer roots of resultant" roots
    if isempty(roots)
        return Poly(one(T), p.var) // 1, Pair(FallingFactorial(p), FallingFactorial(q))
    end

    r = convert(Int64, roots[1])
    g = gcd(shift(p, r), q)
    @debug "GCD" g shift(p, r) q
    if r < 0
        @debug "" [polyval(g, i) for i in 1:-r]
        @debug "" [shift(g, i) for i in 1:-r]
        u = prod(polyval(g, i) for i in 1:-r)
        v = prod(shift(g, i) for i in 1:-r)
        rf = v // u 
    else
        @debug "" [polyval(g, -i+1) for i in 1:r]
        @debug "" [shift(g, -i+1) for i in 1:r]
        u = prod(polyval(g, -i+1) for i in 1:r)
        v = prod(shift(g, -i+1) for i in 1:r)
        rf = u // v
    end
    @debug "" u v rf

    s, sr = divrem(p, shift(g, -r))
    t, tr = divrem(q, g)
    if !iszero(sr) || !iszero(tr)
        @error "Remainder not zero (should not happen) - got $((sr, tr))"
    end
    rfunc, ffact = commonfactors(s, t)
    rfunc * rf, ffact
end

Base.isinteger(x::Basic) = isinteger(convert(Float64, x))

function hgterms(s::RationalFunction)
    c, num, den = monic(s)
    @debug "" c num den
    # fnum, fden = factors(num), factors(den)
    if num != 1 && den != 1
        # TODO: is this shift really the correct thing to do?!
        proots = keys(mroots(num)) |> collect
        qroots = keys(mroots(den)) |> collect
        pqroots = filter(x -> isinteger(x), [proots; qroots])
        @debug "" pqroots
        sh = isempty(pqroots) ? 0 : minimum(pqroots)
        sh = sh < 0 ? -sh : 0
        num, den = shift(num, sh), shift(den, sh)
        @debug "" num den
        rfunc, ffact = commonfactors(num, den)
        @debug "Common factors" rfunc ffact 
        rfunc, ffact = shift(rfunc, -sh), Pair(shift(ffact[1], -sh), shift(ffact[2], -sh))
        @debug "Common factors" rfunc ffact 


        return c, rfunc, ffact
    end
    c, 1, Pair(FallingFactorial(num), FallingFactorial(den))
end

function petkovsek(cf::Vector{T}, arg::T) where {T}
    cf = simplify.(cf)
    @debug "Petkovsek - input" cf
    ls = summands.(cf) |> Iterators.flatten
    ds = Base.denominator.(ls)
    val = lcm2(ds...)
    @info "" ds val
    cf *= val
    @info "" cf
    cf = simplify.(cf)
    cf = coeffs.(cf, arg) # SymPy.Poly.(coeffs, arg)
    hyper = alghyper(Poly.(cf, string(arg)), arg)
    @debug "Petkovsek - alghyper" Base.unique(hyper)
    hgterms.(hyper)
end

# function hyper()

# end # module

# using SymPy
# using Petkovsek

# @syms n y

# @syms n; Recurrences.petkovsek([2*n*(n+1), -(n^2 +3*n-2), n-1], n)

# # algpoly([3, -n, n-1], 0*n, n)
# # algpoly([n-1, -n, 3], 0*n, n)
# println("alghyper: ", [tohg(sol, n) for sol in alghyper([2*n*(n+1), -(n^2 +3*n-2), n-1], n)])
# println("algpoly: ", algpoly([n*(n + 1), -1*n^2 - 3*n + 2, 2*n - 2], 0*n,n))