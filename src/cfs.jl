abstract type ClosedForm end

struct PSolvClosedForm{T} <: ClosedForm
    func::T
    arg::T
    evec::Vector{T}                     # bases of exponentials
    xvec::Vector{T}                     # coeffs
    rvec::Vector{RationalFunction{T}}   # rational functions
    fvec::Vector{RationalFactorial{T}}  # falling factorials
    ivec::Vector{T}                     # initial variables
    instance::T                         # instantiate closed form, yields a closed form where `arg` is replaced by `instance`
end

function CFiniteClosedForm(func::T, arg::T, roots::Vector{T}, coeffs::Vector{T}, ivars::Vector{T}) where {T}
    l = length(roots)
    rs = ones(RationalFunction{T}, l)
    fs = ones(RationalFactorial{T}, l)
    PSolvClosedForm(func, arg, roots, coeffs, rs, fs, ivars, arg)
end

ClosedForm(func::T, c::PSolvClosedForm{T}) where {T} = PSolvClosedForm(func, c.arg, c.evec, c.xvec, c.rvec, c.fvec, c.ivec, c.instance)

func(c::PSolvClosedForm) = c.func
arg(c::PSolvClosedForm) = c.arg
exponentials(c::PSolvClosedForm) = c.evec
coeffs(c::PSolvClosedForm) = c.xvec
rationalfunctions(c::PSolvClosedForm) = c.rvec
factorials(c::PSolvClosedForm) = c.fvec
initvars(c::PSolvClosedForm) = c.ivec

Base.zero(c::PSolvClosedForm{T}) where {T} = PSolvClosedForm(func(c), arg(c), T[], T[], RationalFunction{T}[], RationalFactorial{T}[], T[], arg(c))
Base.iszero(c::PSolvClosedForm) = isempty(exponentials(c)) && isempty(coeffs(c)) && isempty(rationalfunctions(c)) && isempty(factorials(c))

# ------------------------------------------------------------------------------

function Base.:*(c::PSolvClosedForm, x::Number)
    xvec = c.xvec * x
    PSolvClosedForm(c.func, c.arg, c.evec, xvec, c.rvec, c.fvec, c.ivec, c.instance)
end
Base.:*(x::Number, c::ClosedForm) = c * x

function Base.:+(c1::PSolvClosedForm{T}, c2::PSolvClosedForm{T}) where {T}
    @assert arg(c1) == arg(c2) "Argument mismatch, got $(arg(c1)) and $(arg(c1))"
    c1, c2 = reset(c1), reset(c2)
    evec = [c1.evec; c2.evec]
    rvec = [c1.rvec; c2.rvec]
    fvec = [c1.fvec; c2.fvec]
    xvec = [c1.xvec; c2.xvec]
    ivec = Base.unique([c1.ivec; c2.ivec])
    PSolvClosedForm(func(c1), arg(c1), evec, xvec, rvec, fvec, ivec, arg(c1))
end
Base.:-(c1::ClosedForm, c2::ClosedForm) where {T} = c1 + (-1) * c2

# ------------------------------------------------------------------------------

function (c::PSolvClosedForm{T})(n::Union{Int,T}) where {T}
    PSolvClosedForm(c.func, c.arg, c.evec, [subs(x, c.arg, n) for x in c.xvec], c.rvec, c.fvec, c.ivec, subs(c.instance, c.arg, n))
end

init(c::PSolvClosedForm, d::Dict) = PSolvClosedForm(c.func, c.arg, c.evec, [subs(x, d...) for x in c.xvec], c.rvec, c.fvec, c.ivec, c.instance)

function reset(c::PSolvClosedForm{T}) where {T}
    if c.arg in free_symbols(c.instance)
        shift = c.instance - c.arg
        factors = [e^shift for e in c.evec]
        xvec = c.xvec .* factors
        evec = c.evec
    else
        xvec = c.xvec .* (c.evec.^c.instance)
        evec = fill(one(T), length(c.rvec))
    end
    PSolvClosedForm(c.func, c.arg, evec, xvec, c.rvec, c.fvec, c.ivec, c.arg)
end

# ------------------------------------------------------------------------------

function rhs(::Type{T}, c::PSolvClosedForm{T}; expvars = nothing, factvars = nothing) where {T}
    n = c.instance
    exps = expvars
    if expvars == nothing
        exps = exponentials(c).^n
    end
    facts = factvars
    if factvars == nothing
        facts = factorials(c)
    end
    # @info "" c.instance expand.(coeffs(c))
    terms = zip(exps, coeffs(c), rationalfunctions(c), facts)
    # sum(e for (e, x, r, f) in terms) |> expand
    sum(e * x * convert(T, r) * convert(T, f) for (e, x, r, f) in terms) |> expand
end

rhs(::Type{Expr}, c::PSolvClosedForm{T}) where {T} = convert(Expr, rhs(T, c))

function lhs(::Type{Expr}, c::ClosedForm)
    func = Symbol(string(c.func))
    arg = Symbol(string(c.instance))
    :($func($arg))
end

asfunction(c::ClosedForm) = :($(lhs(Expr, c)) = $(rhs(Expr, c)))

# ------------------------------------------------------------------------------

function closedform(rec::CFiniteRecurrence{T}) where {T}

    # TODO: allow inhomogeneous recurrences?
    roots = Poly(coeffs(rec)) |> mroots
    # @info "Roots" roots

    size = order(rec)
    mvec = [T(i) for (_, m) in roots for i in 0:m - 1] # multiplicities
    rvec = [z for (z, m) in roots for _ in 0:m - 1] # roots
    @debug "Roots of characteristic polynomial" collect(zip(rvec, mvec))

    A = [i^m * r^i for i in 0:size - 1, (r, m) in zip(rvec, mvec)]
    b = [initvar(rec.func, i) for i in 0:size - 1] 
    # @info "Ansatz" A b A\b
    # @info "Ansatz" A b linsolve(A, b)
    sol = linsolve(A, b)
    xvec = [x * rec.arg^m for (x, m) in zip(sol, mvec)]
    CFiniteClosedForm(rec.func, rec.arg, rvec, xvec, b)
end

function closedform(rec::HyperRecurrence{T}) where {T}
    hgterms = petkovsek(rec.coeffs, rec.arg)
    
    evec = T[]
    rvec = RationalFunction{T}[]
    fvec = RationalFactorial{T}[]
    for (exp, rfunc, fact) in hgterms
        @debug "" exp rfunc fact
        push!(evec, exp)
        push!(rvec, rfunc)
        push!(fvec, fact)
    end

    size = order(rec)
    A = [e^i * r(i) * f[1](i) / f[2](i) for i in 0:size - 1, (e, r, f) in zip(evec, rvec, fvec)]
    b = [initvar(rec.func, i) for i in 0:size - 1] 
    # @info "" A b linsolve(A, b) hgterms
    PSolvClosedForm(rec.func, rec.arg, evec, linsolve(A, b), rvec, fvec, b, rec.arg)
end

# ------------------------------------------------------------------------------

Base.show(io::IO, c::ClosedForm) = print(io, string(asfunction(c)))

function Base.show(io::IO, ::MIME"text/plain", c::ClosedForm)
    summary(io, c)
    println(io, ":")
    show(io, c)
end