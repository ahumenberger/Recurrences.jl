export Recurrence, LinearRecurrence, CFiniteRecurrence, CFiniteClosedForm
export closedform, expression

abstract type Recurrence end
abstract type LinearRecurrence <: Recurrence end

struct CFiniteRecurrence{T} <: LinearRecurrence
    func::T
    arg::T
    coeffs::Vector{T}
    inhom::T

    function CFiniteRecurrence(func::T, arg::T, coeffs::Vector{T}, inhom::T = T(0)) where {T}
        if any(has.(coeffs, arg))
            error("Not a C-finite recurrence.")
        end
        new{T}(func, arg, coeffs, inhom)
    end
end

struct HyperRecurrence{T} <: LinearRecurrence
    func::T
    arg::T
    coeffs::Vector{T}
    inhom::T

    function HyperRecurrence(func::T, arg::T, coeffs::Vector{T}, inhom::T = T(0)) where {T}
        # if any(has.(coeffs, arg))
        #     error("Not a C-finite recurrence.")
        # end
        new{T}(func, arg, coeffs, inhom)
    end
end

coeffs(r::CFiniteRecurrence) = r.coeffs
order(r::CFiniteRecurrence) = length(r.coeffs) - 1

function Base.show(io::IO, r::CFiniteRecurrence)
    res = join(["$(c) * $(r.func)($(r.arg + (i - 1)))" for (i, c) in enumerate(r.coeffs)], " + ")
    print(io, "$(res) = $(r.inhom)")
end

struct CFiniteClosedForm{T}
    func::T
    arg::T
    mvec::Vector{T}
    rvec::Vector{T}
    xvec::Vector{T}
    initvec::Vector{T}
    instance::T # instantiate closed form, yields a closed form where `arg` is replaced by `instance`
end

CFiniteClosedForm(func::T, arg::T, mvec::Vector{T}, rvec::Vector{T}, xvec::Vector{T}, initvec::Vector{T}) where {T} = CFiniteClosedForm(func, arg, mvec, rvec, xvec, initvec, arg)
CFiniteClosedForm(func::T, cf::CFiniteClosedForm{T}) where {T} = CFiniteClosedForm(func, cf.arg, cf.mvec, cf.rvec, cf.xvec, cf.initvec, cf.instance)

function Base.:*(cf::CFiniteClosedForm{T}, coeff::Number) where {T}
    xvec = coeff .* cf.xvec
    CFiniteClosedForm(cf.func, cf.arg, cf.mvec, cf.rvec, xvec, cf.initvec, cf.instance)
end
Base.:*(coeff::Number, cf::CFiniteClosedForm{T}) where {T} = cf*coeff

function Base.:+(cf1::CFiniteClosedForm{T}, cf2::CFiniteClosedForm{T}) where {T}
    @assert cf1.arg == cf2.arg "Argument mismatch, got $(cf1.arg) and $(cf2.arg)"
    cf1, cf2 = reset(cf1), reset(cf2)
    mvec = [cf1.mvec; cf2.mvec]
    rvec = [cf1.rvec; cf2.rvec]
    xvec = [cf1.xvec; cf2.xvec]
    initvec = [cf1.initvec; cf2.initvec]
    CFiniteClosedForm(cf1.func, cf1.arg, mvec, rvec, xvec, initvec)
end
Base.:-(cf1::CFiniteClosedForm{T}, cf2::CFiniteClosedForm{T}) where {T} = cf1 + (-1) * cf2

function (c::CFiniteClosedForm{T})(n::Union{Int, T}) where {T}
    CFiniteClosedForm(c.func, c.arg, c.mvec, c.rvec, [subs(x, c.arg, n) for x in c.xvec], c.initvec, subs(c.instance, c.arg, n))
end

function reset(c::CFiniteClosedForm{T}) where {T}
    if has(c.instance, c.arg)
        shift = c.instance - c.arg
        factors = [c.instance^m / c.arg^m * r^shift for (r, m) in zip(c.rvec, c.mvec)]
        xvec = c.xvec .* factors
        rvec = c.rvec
    else
        xvec = c.xvec .* (c.rvec .^ c.instance)
        rvec = fill(T(1), length(c.rvec))
    end
    CFiniteClosedForm(c.func, c.arg, c.mvec, rvec, xvec, c.initvec)
end

init(c::CFiniteClosedForm, d::Dict) = CFiniteClosedForm(c.func, c.arg, c.mvec, c.rvec, [subs(x, d) for x in c.xvec], c.initvec, c.instance)

function mroots(poly::Poly{T}) where {T}
    roots = Polynomials.roots(poly)
    Dict([(r, count(x -> x==r, roots)) for r in Base.unique(roots)])
end

function mroots(p::Poly{SymPy.Sym})
    T = SymPy.Sym
    length(p) == 0 && return zeros(T, 0)

    num_leading_zeros = 0
    while p[num_leading_zeros] == zero(T)
        if num_leading_zeros == length(p)-1
            return zeros(T, 0)
        end
        num_leading_zeros += 1
    end
    num_trailing_zeros = 0
    while p[end - num_trailing_zeros] == zero(T)
        num_trailing_zeros += 1
    end
    n = lastindex(p)-(num_leading_zeros + num_trailing_zeros)
    n < 1 && return zeros(T, length(p) - num_trailing_zeros - 1)

    companion = diagm(-1 => ones(T, n-1))
    an = p[end-num_trailing_zeros]
    companion[1,:] = -p[(end-num_trailing_zeros-1):-1:num_leading_zeros] / an

    companion[:eigenvals]()
end

function closedform(rec::CFiniteRecurrence{T}) where {T}

    # TODO: allow inhomogeneous recurrences?
    roots = Poly(coeffs(rec)) |> mroots
    # @info "Roots" roots

    size = order(rec)
    mvec = [T(i) for (_, m) in roots for i in 0:m - 1] # multiplicities
    rvec = [z for (z, m) in roots for _ in 0:m - 1] # roots
    @debug "Roots of characteristic polynomial" collect(zip(rvec, mvec))

    A = [i^m * r^i for i in 0:size-1, (r, m) in zip(rvec, mvec)]
    b = [initvariable(rec.func, i) for i in 0:size - 1] 
    # @info "Ansatz" A b A\b
    CFiniteClosedForm(rec.func, rec.arg, mvec, rvec, A \ b, b)
end

function expression(cf::CFiniteClosedForm)
    vec = [cf.instance^m * r^cf.instance for (r, m) in zip(cf.rvec, cf.mvec)]
    simplify(transpose(vec) * cf.xvec)
end

Base.show(io::IO, cf::CFiniteClosedForm) = print(io, " $(cf.func)($(cf.instance)) = $(expression(cf))")

function Base.show(io::IO, ::MIME"text/plain", cf::CFiniteClosedForm)
    summary(io, cf)
    println(io, ":")
    show(io, cf)
end

struct HyperClosedForm{T}
    func::T
    arg::T
    evec::Vector{T} # bases of exponentials
    rvec::Vector{RationalFunction{T}} # rational functions
    fvec::Vector{Pair{FallingFactorial{T},FallingFactorial{T}}} # falling factorials
    initvec::Vector{T}
    instance::T # instantiate closed form, yields a closed form where `arg` is replaced by `instance`
end


function closedform(rec::HyperRecurrence{T}) where {T}
    hgterms = petkovsek(rec.coeffs, rec.arg)
    
    evec = T[]
    rvec = RationalFunction{T}[]
    fvec = Vector{Pair{FallingFactorial{T},FallingFactorial{T}}}[]
    for (exp, rfunc, fact) in hgterms
        push!(evec, exp)
        push!(rvec, rfunc)
        push!(fvec, fact)
    end

    # A = [e^i * r(i) * f[1](i) / f[2](i) for i in 0:size-1, (e, r, f) in zip(evec, rvec, fvec)]
    # b = [initvariable(rec.func, i) for i in 0:size - 1] 

    # HyperClosedForm(rec.func, rec.arg, evec, rvec, fvec)
end