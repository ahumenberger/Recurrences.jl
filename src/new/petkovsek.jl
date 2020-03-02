function algpoly(plist::Vector{T}, f, n) where {T}
    @debug "Algorithm [algpoly]" plist f
    R = parent(n)
    # Construct polynomials qj that arise from treating the shift operator
	# as the difference operator plus the identity operator.
    r = length(plist) - 1
    qlist = [sum(binomial(i, j) * plist[i+1] for i in j:r) for j in 0:r]

    # Find all candidates for the degree bound on output polynomial.
    b = maximum([(Nemo.degree(q) - (j - 1)) for (j, q) in enumerate(qlist)])
    alpha = sum(lead(q)*fallfactorial(n, j-1) for (j, q) in enumerate(qlist) if Nemo.degree(q) - (j-1) == b)
    roots = [-Nemo.coeff(f, 0) for (f, _) in Nemo.factor(R(alpha))]
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

function monic_factors(p::PolyElem)
    fac = Nemo.factor(p)
    fs = [(lead(f), div(f, parent(p)(lead(f)))) for (f, m) in fac for i in 1:m]
    Nemo.coeff(parent(p)(unit(fac)), 0)*prod(map(first, fs)), map(last, fs)
end

function all_monic_factors(p::PolyElem)
    unit, fs = monic_factors(p)
    fs = [one(p); fs]
    unique(prod(s) for s in Combinatorics.powerset(fs, 1))
end

_prod(itr) = reduce(*, itr, init=1)

function alghyper(plist::Vector{T}, n; all_solutions::Bool=false) where {T}
    @debug "Algorithm [alghyper]" plist
    d = length(plist)

    alist = all_monic_factors(plist[1](n))
    blist = all_monic_factors(plist[end](n-(d-2)))

    @debug "All monic factors" alist blist

    solutions = FracElem{typeof(n)}[]
    for a in alist
        for b in blist
            @debug "Monic factors" a b

            ps = [plist[i+1] * _prod(a(n+j) for j in 0:i-1) * _prod(b(n+j) for j in i:d-2) for i in 0:d-1]

            m = maximum(Nemo.degree(p) for p in ps)
            alphas = map(x->Nemo.coeff(x, m), ps)
            zpoly = sum(c*n^(i-1) for (i, c) in enumerate(alphas))
            zs = roots(zpoly)
            zs = [z for z in zs if !iszero(z)]
            @debug "Nonzero roots" alphas zpoly zs
            for z in zs
                _ps = [z^(i-1)*p for (i, p) in enumerate(ps)]
                @debug "Coefficients auxiliary recurrence" _ps
                
                _sols = algpoly(_ps, zero(n), n)
                @debug "Solutions of algpoly" _sols

                for c in _sols
                    s = z * (a(n) // b(n)) * (c(n+1) // c(n))
                    !all_solutions && return s
                    push!(solutions, s)
                end
            end
        end
    end
    @debug "Algorithm [alghyper] return" unique(solutions) solutions
    !all_solutions && return nothing
    return unique(solutions)
end

function gosper(t::Seq{T}, arg) where {T <: FieldElem}
    @debug "[Gosper]" t arg
    r = div(t(arg+1), t(arg))
    @assert isrational(r)
    r = first(coeffs(r))
    c1, fac1 = monic_factors(numerator(r))
    c2, fac2 = monic_factors(denominator(r))
    if isempty(fac1) && isempty(fac2)
        @debug "Rational coefficients" c1 c2
        return parent(t)(c1//c2*arg)
    end
    R = parent(arg)
    f = R(_prod(fac1))
    g = R(_prod(fac2))
    @debug "Gosper start" t r f g
    @assert gcd(f, g) == 1

    _R, _h = PolynomialRing(base_ring(R), "h")
    _arg = change_base_ring(_R, arg)
    _res = resultant(f(_arg), g(_arg+_h))
    @debug "Resultant" _res

    hs = roots(_res)
    hs = [numerator(r) for r in hs if r >= 0 && isone(denominator(r))]
    @debug "Non-negative integer roots of resultant" hs

    sz, type = length(hs), typeof(arg)
    ss = Array{type}(undef, sz)
    ps = Array{type}(undef, sz+1)
    qs = Array{type}(undef, sz+1)
    ps[1] = f
    qs[1] = g
    for j in 1:length(hs)
        ss[j] = gcd(ps[j](arg), qs[j](arg+hs[j]))
        ps[j+1] = div(ps[j], ss[j])
        qs[j+1] = div(qs[j], ss[j](arg-hs[j]))
    end
    a = c1//c2 * ps[end]
    b = qs[end]
    c = R(_prod(ss[i](arg-j) for i in 1:length(hs) for j in 1:hs[i]))
    @debug "Coefficients for auxiliary equation" a b c ss

    if degree(a) != degree(b) || lead(a) != lead(b)
        D = [degree(c) - max(degree(a), degree(b))]
    else
        d = degree(a)
        if iszero(d)
            D = [degree(c) - d + 1]
        else
            A = Nemo.coeff(a, d-1)
            B = Nemo.coeff(b(arg-1), d-1)
            dd = (B-A)/lead(a)
            if dd >= 0 && isone(denominator(dd))
                D = [degree(c) - d + 1, convert(Int, numerator(dd))]
            else
                D = [degree(c) - d + 1]
            end
        end
    end
    filter!(x->x>=0, D)
    @debug "Degree set" D
    isempty(D) && return nothing
    d = maximum(D)

    # Use method of undetermined coefficients
    _S, _vs = PolynomialRing(base_ring(R), ["v$i" for i = 0:d])
    _arg = change_base_ring(_S, arg)
    _a = change_base_ring(_S, a)
    _b = change_base_ring(_S, b)
    aux = sum(_vs[i+1]*_arg^i for i = 0:d)
    poly = _a(_arg)*aux(_arg+1) - _b(_arg-1)*aux(_arg)
    @debug "Solution template" aux poly

    deg = degree(poly)
    cs = [Nemo.coeff(poly, i) for i in 0:deg]
    mA = [Nemo.coeff(c, v) for c in cs, v in _vs]
    mB = [Nemo.coeff(c, i) for i in 0:deg]
    @debug "Solve" mA mB

    mM = hcat(mA, mB)
    mM = MatrixSpace(base_ring(arg), size(mM)...)(mM)
    r, mM = rref(mM)

    xs = zeros(base_ring(arg), d+1)
    nc = ncols(mM)
    for i in 1:r
        j = findfirst(isone, collect(mM[i, 1:nc-1]))
        xs[j] = mM[i, nc]
    end
    @debug "Method of undetermined coefficients" r mM xs
    iszero(xs) && return nothing

    x = sum(xs[i+1]*arg^i for i = 0:d)
    z = parent(t)(b(arg-1)*x(arg)//c(arg))*t
    @debug "Solution from gosper" z x b c
    @assert z(arg+1) - z(arg) == t
    return z
end

issolution(s::Seq, cfs::Vector{<:FracElem}, arg::PolyElem) = iszero(sum(c * s(arg + i - 1) for (i, c) in enumerate(cfs)))

function petkovsek(cfs::Vector{<:FracElem}, arg::PolyElem, level::Int=0)
    @debug "[Petkovsek]" cfs arg
    # Get rid of denominator
    m = reduce(lcm, map(denominator, cfs))
    polys = [numerator(c * m) for c in cfs]

    ls = map(monic_factors, polys)
    cs = map(first, ls)
    g = reduce(gcd, cs)
    cs = map(x->div(x, g), cs)
    polys = [parent(arg)(c*p) for (c, p) in zip(cs, map(_prod∘last, ls))]


    g = reduce(gcd, polys)
    polys = map(x->div(x, g), polys)
    S = SequenceRing(parent(arg))
    sz = length(polys)
    if sz == 2
        # Recurrence is hypergeometric
        return [S((-polys[1]//polys[2])(arg), arg+1)(arg-1)]
    end

    # Get first hypergeometric solution
    hgterm = alghyper(polys, arg)
    if isnothing(hgterm)
        @debug "No solution from Hyper"
        return nothing
    end
    s_n = S(hgterm(arg), arg+1)(arg-1)
    @debug "Hyper" hgterm s_n
    @assert issolution(s_n, cfs, arg)
    
    # Take y(n) = z(n)s(n) for unknown z(n), and plug into original recurrence.
    # Then reduce order by taking u(n) = z(n+1) - z(n)
    R, zs = PolynomialRing(S, ["z$i" for i in 0:sz-1])
    poly = sum(s_n(arg+(i-1))*S(cfs[i])*zs[i] for i in 1:sz)

    bs = [Nemo.coeff(poly, z) for z in zs]
    for i in reverse(2:length(bs)-1)
        bs[i] = bs[i] + bs[i+1]
    end
    bs = bs[2:end]
    g = reduce(gcd, bs)
    bs = map(x->div(x, g), bs)
    @debug "Coefficients for order-reduced recurrence" g bs map(isrational, bs) terms(bs[end])
    @assert all(isrational(b) for b in bs)
    # Recursively solve for u(n)
    u_basis = petkovsek(map(first∘coeffs, bs), arg, level+1)
    @assert all(iszero(sum(c * u(arg + i - 1) for (i, c) in enumerate(bs))) for u in u_basis)
    # Use Gosper to resolve antidifference, i.e. find z(n) satisfying u(n) = z(n+1) - z(n)
    z_basis = [gosper(u, arg) for u in u_basis]
    if any(map(isnothing, z_basis))
        @debug "No solution from Gosper"
        return nothing
    end
    @debug "Solutions from Gosper" z_basis
    ss = [s_n * z for z in z_basis]
    @assert all(issolution(s, cfs, arg) for s in ss)
    [s_n; ss]
end