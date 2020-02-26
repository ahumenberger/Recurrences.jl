using Nemo
using Hecke
using Combinatorics

fallfactorial(n, j) = j == 0 ? 1 : prod(n - i for i in 0:j-1)

function algpoly(plist::Vector{T}, f, n) where {T}
    @debug "Algorithm [algpoly]" plist f
    # Construct polynomials qj that arise from treating the shift operator
	# as the difference operator plus the identity operator.
    r = length(plist) - 1
    qlist = [sum(binomial(i, j) * plist[i+1] for i in j:r) for j in 0:r]

    # Find all candidates for the degree bound on output polynomial.
    b = maximum([(Nemo.degree(q) - (j - 1)) for (j, q) in enumerate(qlist)])
    alpha = sum(lc(q)*fallfactorial(n, j-1) for (j, q) in enumerate(qlist) if Nemo.degree(q) - (j-1) == b)
    roots = [-Nemo.coeff(f, 0) for (f, _) in Nemo.factor(alpha)]
    roots = [r for r in roots if r >= 0 && isone(denominator(r))]

    first = Nemo.degree(f) - b
    second = -b - 1
    third = maximum([roots; 0])
    d = convert(Int, numerator(max(first, second, third)))
    @debug "Degree bound" first second third

    # Use method of undetermined coefficients to find the output polynomial.
    R, varlist = PolynomialRing(base_ring(n), ["x$i" for i in 1:d+1])
    arg = change_base_ring(R, n)
    genpoly = sum(v*arg^(i-1) for (i,v) in enumerate(varlist))
    _plist = map(x->change_base_ring(R, x), plist)
    @debug "Generic polynomial" genpoly varlist

    poly = sum(genpoly(arg+(i - 1)) * p for (i, p) in enumerate(_plist))
    poly -= change_base_ring(R, f)
    @debug "Poly" poly

    coefficients = [Nemo.coeff(poly, i) for i in 0:length(poly)]
    A = [Nemo.coeff(cf, v) for cf in coefficients, v in varlist]
    A = MatrixSpace(base_ring(R), size(A)...)(A)
    rank, ker = kernel(A)
    solutions = typeof(n)[]
    S = base_ring(n)
    for j in 1:rank
        col = ker[:, j]
        s = sum(Nemo.coeff(genpoly, i)(col...)*n^i for i in 0:length(genpoly))
        push!(solutions, s)
    end
    @debug "Algorithm [algpoly] return" solutions
    return solutions
end

function factors(p::fmpq_poly)
    fac = Nemo.factor(p)
    unit(fac), first.(collect(fac))
end

factor_list(p::fmpq_poly) = factors(p)[2]
lc(p::fmpq_poly) = factors(p)[1]

function monic_factors(p::PolyElem)
    fs = [f for (f, m) in Nemo.factor(p) for i in 1:m]
    fs = [one(p); fs]
    unique(prod(s) for s in Combinatorics.powerset(fs, 1))
end

function _roots(p::PolyElem)
    [-Nemo.coeff(f, 0) for (f, _) in Nemo.factor(p)]
end

_prod(itr) = reduce(*, itr, init=1)

function alghyper(plist::Vector{T}, n) where {T}
    @debug "Algorithm [alghyper]" plist
    d = length(plist)

    alist = monic_factors(plist[1](n))
    blist = monic_factors(plist[end](n-(d-2)))

    @debug "All monic factors" alist blist

    solutions = FracElem{typeof(n)}[]
    for a in alist
        for b in blist
            @debug "Monic factors" a b

            ps = [plist[i+1] * _prod(a(n+j) for j in 0:i-1) * _prod(b(n+j) for j in i:d-1) for i in 0:d-1]

            m = maximum(Nemo.degree(p) for p in ps)
            alphas = map(x->Nemo.coeff(x, m), ps)
            zpoly = sum(c*n^(i-1) for (i, c) in enumerate(alphas))
            zs = _roots(zpoly)
            zs = [z for z in zs if !iszero(z)]
            @debug "Nonzero roots" alphas zpoly zs
            for z in zs
                _ps = [z^(i-1)*p for (i, p) in enumerate(ps)]
                @debug "Coefficients auxiliary recurrence" _ps
                
                _sols = algpoly(_ps, zero(n), n)
                @debug "Solutions of algpoly" _sols

                for c in _sols
                    s = z * (a(n) // b(n)) * (c(n+1) // c(n))
                    push!(solutions, s)
                end
            end
        end
    end
    @debug "Algorithm [alghyper] return" unique(solutions)
    return unique(solutions)
end

# function commonfactors(p::Poly{T}, q::Poly{T}) where {T}
#     if p == 1 || q == 1
#         return Poly(one(T), p.var) // 1, Pair(FallingFactorial(p), FallingFactorial(q))
#     end
#     k = variables(T)
#     res = resultant(shift(p, k), q)
#     roots = keys(mroots(res)) |> collect
#     filter!(x -> isinteger(x), roots)
#     @debug "Integer roots of resultant" roots
#     if isempty(roots)
#         return Poly(one(T), p.var) // 1, Pair(FallingFactorial(p), FallingFactorial(q))
#     end

#     r = convert(Int64, roots[1])
#     g = gcd(shift(p, r), q)
#     @debug "GCD" g shift(p, r) q
#     if r < 0
#         @debug "" [polyval(g, i) for i in 1:-r]
#         @debug "" [shift(g, i) for i in 1:-r]
#         u = prod(polyval(g, i) for i in 1:-r)
#         v = prod(shift(g, i) for i in 1:-r)
#         rf = v // u 
#     else
#         @debug "" [polyval(g, -i+1) for i in 1:r]
#         @debug "" [shift(g, -i+1) for i in 1:r]
#         u = prod(polyval(g, -i+1) for i in 1:r)
#         v = prod(shift(g, -i+1) for i in 1:r)
#         rf = u // v
#     end
#     @debug "" u v rf

#     s, sr = divrem(p, shift(g, -r))
#     t, tr = divrem(q, g)
#     if !iszero(sr) || !iszero(tr)
#         @error "Remainder not zero (should not happen) - got $((sr, tr))"
#     end
#     rfunc, ffact = commonfactors(s, t)
#     rfunc * rf, ffact
# end

# function hgterms(s::RationalFunction)
#     c, num, den = monic(s)
#     @debug "" c num den
#     # fnum, fden = factors(num), factors(den)
#     if num != 1 && den != 1
#         # TODO: is this shift really the correct thing to do?!
#         proots = keys(mroots(num)) |> collect
#         qroots = keys(mroots(den)) |> collect
#         pqroots = filter(x -> isinteger(x), [proots; qroots])
#         @debug "" pqroots
#         sh = isempty(pqroots) ? 0 : minimum(pqroots)
#         sh = sh < 0 ? -sh : 0
#         num, den = shift(num, sh), shift(den, sh)

#         rfunc, ffact = commonfactors(num, den)
#         rfunc, ffact = shift(rfunc, -sh), Pair(shift(ffact[1], -sh), shift(ffact[2], -sh))

#         return c, rfunc, ffact
#     end
#     c, 1, Pair(FallingFactorial(num), FallingFactorial(den))
# end

function petkovsek(coeffs::Vector{T}, arg::S) where {S,T}
    @info "Petkovsek - input" coeffs
    ds = map(denominator, coeffs)
    p = reduce(lcm, ds)
    polys = reverse([numerator(c * p) for c in coeffs])
    @info "" polys ds p
    hyper = alghyper(polys, arg)
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