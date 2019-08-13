using Nemo
using Hecke

function algpoly(polylist::Vector{T}, f, n) where {T}
    @info "algpoly" polylist
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
    b = maximum([(Nemo.degree(poly) - (j - 1)) for (j, poly) in enumerate(qlist)])
    first = Nemo.degree(f) - b
    second = -1 * b - 1
    lcoeffs = lc.(qlist)
    alpha = 0
    for j in 0:length(qlist)-1
        aux = Nemo.degree(qlist[j+1]) - j
        if (Nemo.degree(qlist[j+1]) - j) == b
            alpha += lcoeffs[j+1] * fallingfactorial(n, j)
        end
    end
    @debug "" alpha typeof(alpha)
    deg = Nemo.roots(alpha)
    deg = isa(deg, Dict) ? keys(deg) : deg
    third = maximum(deg)
    @info "" qlist b deg first second third
    d = convert(Int64, denominator(max(first, second, third, 0)))

    # Use method of undetermined coefficients to find the output polynomial.
    varlist = mkvar.(["x$i" for i in 1:d+1])
    arg = mkvar(string(n))
    # pol = 0 * n
    # for i in 0:d
    #     pol += n^i * varlist[i+1]
    # end
    # R, varlist = Nemo.PolynomialRing(Nemo.QQ, [vars; string(n)])
    @info "" varlist
    genpoly = sum(v*arg^(i-1) for (i,v) in enumerate(varlist))
    @info "Generic polynomial" genpoly d+1 varlist
    # solution = -1 * f
    # for (i, p) in enumerate(polylist)
    #     @debug "" shift(poly, i - 1) * p
    #     # solution += subs(pol, (n, n + (i-1))) * poly
    #     # @debug "" solution
    # end

    poly = sum(MultivariatePolynomials.subs(genpoly, (arg=>arg+(i - 1))) * mkpoly(string(p)) for (i, p) in enumerate(polylist))
    @info "" poly  [MultivariatePolynomials.subs(genpoly, (arg=>arg+(i - 1))) * mkpoly(string(p)) for (i, p) in enumerate(polylist)]
    poly -= mkpoly(string(f))
    @info "Poly" poly
    if iszero(poly)
        return [parent(f)(1)]
    end
    coefficients = coeffs(poly)
    filter!(e -> e != 0, coefficients)
    sol = solve(coefficients, varlist)
    missing = setdiff(varlist, keys(sol))
    p = Poly([subs(c, sol) for c in coeffs(genpoly)], string(n))
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

function factors(p::fmpq_poly)
    fac = Nemo.factor(p)
    unit(fac), first.(collect(fac))
end

factor_list(p::fmpq_poly) = factors(p)[2]
lc(p::fmpq_poly) = factors(p)[1]

function alghyper(polylist::Vector{T}, n) where {T}
    @info "-> alghyper" polylist n
    p = polylist[1]
    lcoeff, alist = factors(p)
    @info "" alist lcoeff

    # for (i,p) in enumerate(alist)
    #     alist[i] /= lcoeff
    # end
    if !(1 in alist)
        push!(alist, n)
    end

    d = length(polylist)
    p = polylist[end](n+(-d+2))
    @info "" p
    lcoeff = lc(p)
    blist = factor_list(p)
    # for (i,p) in enumerate(blist)
    #     blist[i] /= lcoeff
    # end
    if !(1 in blist)
        push!(blist, n)
    end

    alist, blist = reverse(alist), reverse(blist)

    @info "" alist blist

    solutions = []
    for aelem in alist
        for belem in blist
            @info "Monic factors" aelem belem
            plist = T[]
            for i in 0:d-1 
                pi = polylist[i+1]
                for j in 0:i-1
                    pi *= aelem(n+j)
                end
                for j in i:d-1
                    pi *= belem(n+j)
                end
                push!(plist, pi)
            end

            m = maximum(Nemo.degree.(plist))
            alpha = Nemo.coeff.(plist, m)
            @info "" alpha
            # @syms z
            # zpol = 0*z
            # for i in 0:length(alpha) - 1
            #     zpol += alpha[i+1]*z^i
            # end
            zpol = sum(c*n^(i-1) for (i,c) in enumerate(alpha))

            @info "" zpol typeof(zpol) typeof(alpha) typeof(plist)
            vals = Nemo.roots(zpol) #[key for (key,val) in mroots(zpol) if key != 0]
            @info "Roots" zpol vals
            for x in vals
                polylist2 = [x^(i-1)*p for (i,p) in enumerate(plist)]
                @info "" polylist2
                
                polysols = algpoly(polylist2, 0*n, n)
                @debug "Solutions of algpoly" polysols

                arg = mkvar(string(n))
                for c in polysols
                    @debug "" c
                    fn = mkpoly ∘ string
                    s = fn(x) * (fn(aelem) / fn(belem)) * (MultivariatePolynomials.subs(fn(c), arg=>arg+1) / fn(c))
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

        rfunc, ffact = commonfactors(num, den)
        rfunc, ffact = shift(rfunc, -sh), Pair(shift(ffact[1], -sh), shift(ffact[2], -sh))

        return c, rfunc, ffact
    end
    c, 1, Pair(FallingFactorial(num), FallingFactorial(den))
end

function petkovsek(cf::Vector{T}, arg::S) where {S,T}
    # cf = simplify.(cf)
    @info "Petkovsek - input" cf
    # ls = summands.(cf) |> Iterators.flatten
    # ds = Base.denominator.(ls)
    # val = lcm2(ds...)
    ds = map(denominator, cf)
    val = reduce(lcm, ds)
    @info "" ds val
    cf = RPoly[c * val for c in cf]
    @info "" cf
    R, n = Nemo.PolynomialRing(Nemo.QQ, string(arg))
    nemopolys = map(R ∘ string, cf)
    @info "nemopolys" nemopolys
    # cf = simplify.(cf)
    # cf = map(coefficients, cf)
    # alghyper([R("n^2-1")], arg)
    hyper = alghyper(reverse(nemopolys), n)
    @info "Petkovsek - alghyper" hyper
    hgterms.(hyper)
end

(R::FmpqPolyRing)(s::String) = R(Meta.parse(s))

function (R::FmpqPolyRing)(p::Expr)
    v = gen(R)
    vs = [:($(Symbol(string(v))) = $v)]
    q = quote
        let $(vs...)
            $(p)
        end
    end
    eval(q)
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