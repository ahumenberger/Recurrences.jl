abstract type ClosedForm end

struct PSolvClosedForm{T} <: ClosedForm
    func::T
    arg::T
    evec::Vector{T} # bases of exponentials
    xvec::Vector{T} # coeffs
    rvec::Vector{RationalFunction{T}} # rational functions
    fvec::Vector{Pair{FallingFactorial{T},FallingFactorial{T}}} # falling factorials
    # initvec::Vector{T}
    instance::T # instantiate closed form, yields a closed form where `arg` is replaced by `instance`
end

exponentials(c::PSolvClosedForm) = c.evec
coeffs(c::PSolvClosedForm) = c.xvec
rationalfunctions(c::PSolvClosedForm) = c.rvec
factorials(c::PSolvClosedForm) = c.fvec

Base.zero(c::PSolvClosedForm) = PSolvClosedForm(c.func, c.arg, [], [], [], [], c.arg)
Base.iszero(c::PSolvClosedForm) = isempty(exponentials(c)) && isempty(coeffs(c)) && isempty(rationalfunctions(c)) && isempty(factorials(c))

# ------------------------------------------------------------------------------

import Base: zero


struct CFiniteClosedForm{T} <: ClosedForm
    func::T
    arg::T
    mvec::Vector{T}
    rvec::Vector{T}
    xvec::Vector{T}
    initvec::Vector{T}
    instance::T # instantiate closed form, yields a closed form where `arg` is replaced by `instance`
end

CFiniteClosedForm(func::T, arg::T, mvec::Vector{T}, rvec::Vector{T}, xvec::Vector{T}, initvec::Vector{T}) where {T} = CFiniteClosedForm(func, arg, mvec, rvec, xvec, initvec, arg)

struct HyperClosedForm{T} <: ClosedForm
    func::T
    arg::T
    evec::Vector{T} # bases of exponentials
    rvec::Vector{RationalFunction{T}} # rational functions
    fvec::Vector{Pair{FallingFactorial{T},FallingFactorial{T}}} # falling factorials
    xvec::Vector{T} # coeffs
    initvec::Vector{T}
    instance::T # instantiate closed form, yields a closed form where `arg` is replaced by `instance`
end

HyperClosedForm(func::T, arg::T, evec::Vector{T}, rvec::Vector{RationalFunction{T}}, fvec::Vector{Pair{FallingFactorial{T},FallingFactorial{T}}}, xvec::Vector{T}, initvec::Vector{T}) where {T} = HyperClosedForm(func, arg, evec, rvec, fvec, xvec, initvec, arg)

# ------------------------------------------------------------------------------

ClosedForm(func::T, cf::CFiniteClosedForm{T}) where {T} = CFiniteClosedForm(func, cf.arg, cf.mvec, cf.rvec, cf.xvec, cf.initvec, cf.instance)
ClosedForm(func::T, cf::HyperClosedForm{T}) where {T} = HyperClosedForm(func, cf.arg, cf.evec, cf.rvec, cf.fvec, cf.xvec, cf.initvec, cf.instance)

zero(c::CFiniteClosedForm{T}) where {T} = c * 0
zero(c::HyperClosedForm{T}) where {T} = c * 0

# ------------------------------------------------------------------------------

function closedform(rec::CFiniteRecurrence{T}) where {T}

    # TODO: allow inhomogeneous recurrences?
    roots = Poly(coeffs(rec)) |> mroots
    # @info "Roots" roots

    size = order(rec)
    mvec = [T(i) for (_, m) in roots for i in 0:m - 1] # multiplicities
    rvec = [z for (z, m) in roots for _ in 0:m - 1] # roots
    @debug "Roots of characteristic polynomial" collect(zip(rvec, mvec))

    A = [i^m * r^i for i in 0:size-1, (r, m) in zip(rvec, mvec)]
    b = [initvar(rec.func, i) for i in 0:size - 1] 
    # @info "Ansatz" A b A\b
    CFiniteClosedForm(rec.func, rec.arg, mvec, rvec, A \ b, b)
end

function closedform(rec::HyperRecurrence{T}) where {T}
    hgterms = petkovsek(rec.coeffs, rec.arg)
    
    evec = T[]
    rvec = RationalFunction{T}[]
    fvec = Pair{FallingFactorial{T},FallingFactorial{T}}[]
    for (exp, rfunc, fact) in hgterms
        @debug "" exp rfunc fact
        push!(evec, exp)
        push!(rvec, rfunc)
        push!(fvec, fact)
    end

    size = order(rec)
    A = [e^i * r(i) * f[1](i) / f[2](i) for i in 0:size-1, (e, r, f) in zip(evec, rvec, fvec)]
    b = [initvar(rec.func, i) for i in 0:size - 1] 
    @info "" A b
    HyperClosedForm(rec.func, rec.arg, evec, rvec, fvec, A \ b, b)
end

# ------------------------------------------------------------------------------

exponentials(c::CFiniteClosedForm) = c.rvec

function expression(cf::CFiniteClosedForm; expvars = nothing)
    if expvars == nothing
        vec = [cf.instance^m * r^cf.instance for (r, m) in zip(cf.rvec, cf.mvec)]
    else
        @assert length(expvars) == length(cf.rvec) "Number of variables must be equal to number of exponentials"
        vec = [cf.instance^m * r for (r, m) in zip(expvars, cf.mvec)]
    end
    simplify(transpose(vec) * cf.xvec)
end

function convert(::Type{T}, c::CFiniteClosedForm) where {T <: Union{SymPy.Sym,SymEngine.Basic}}
    vec = [c.instance^m * r^c.instance for (r, m) in zip(c.rvec, c.mvec)]
    simplify(transpose(vec) * c.xvec)
end

function convert(::Type{T}, c::HyperClosedForm{T}) where {T <: Union{SymPy.Sym,SymEngine.Basic}}
    vec = [e^c.arg * convert(T, r) * convert(T, f[1]) / convert(T, f[2]) * x for (e, r, f, x) in zip(c.evec, c.rvec, c.fvec, c.xvec)]
    simplify(sum(vec))
end

function convert(::Type{Expr}, c::ClosedForm)
    rhs = convert(Expr, convert(Basic, c))
    func = Symbol(string(c.func))
    arg = Symbol(string(c.instance))
    :($func($arg) = $rhs)
end

rhs(::Type{Expr}, c::ClosedForm) = convert(Expr, convert(Basic, c))

function lhs(::Type{Expr}, c::ClosedForm)
    func = Symbol(string(c.func))
    arg = Symbol(string(c.instance))
    :($func($arg))
end

asfunction(::Type{Expr}, c::ClosedForm) = :($(lhs(Expr, c)) = rhs(Expr, c))

# function expression(cf::HyperClosedForm)
#     vec = [e^i * r(i) * f[1](i) / f[2](i) for i in 0:size-1, (e, r, f) in zip(evec, rvec, fvec)]
#     vec = [cf.instance^m * r^cf.instance for (r, m) in zip(cf.rvec, cf.mvec)]
#     simplify(transpose(vec) * cf.xvec)
# end

# ------------------------------------------------------------------------------

function Base.:*(cf::CFiniteClosedForm{T}, coeff::Number) where {T}
    xvec = coeff .* cf.xvec
    CFiniteClosedForm(cf.func, cf.arg, cf.mvec, cf.rvec, xvec, cf.initvec, cf.instance)
end
Base.:*(coeff::Number, cf::CFiniteClosedForm{T}) where {T} = cf * coeff

function Base.:+(cf1::CFiniteClosedForm{T}, cf2::CFiniteClosedForm{T}) where {T}
    @assert cf1.arg == cf2.arg "Argument mismatch, got $(cf1.arg) and $(cf2.arg)"
    cf1, cf2 = reset(cf1), reset(cf2)
    mvec = [cf1.mvec; cf2.mvec]
    rvec = [cf1.rvec; cf2.rvec]
    xvec = [cf1.xvec; cf2.xvec]
    initvec = [cf1.initvec; cf2.initvec]
    CFiniteClosedForm(cf1.func, cf1.arg, mvec, rvec, xvec, initvec)
end
Base.:-(cf1::ClosedForm, cf2::ClosedForm) where {T} = cf1 + (-1) * cf2

function Base.:*(c::HyperClosedForm{T}, coeff::Number) where {T}
    xvec = coeff .* c.xvec
    HyperClosedForm(c.func, c.arg, c.evec, c.rvec, c.fvec, xvec, c.initvec, c.instance)
end
Base.:*(coeff::Number, c::HyperClosedForm{T}) where {T} = c * coeff

function Base.:+(cf1::HyperClosedForm{T}, cf2::HyperClosedForm{T}) where {T}
    @assert cf1.arg == cf2.arg "Argument mismatch, got $(cf1.arg) and $(cf2.arg)"
    cf1, cf2 = reset(cf1), reset(cf2)
    evec = [cf1.evec; cf2.evec]
    rvec = [cf1.rvec; cf2.rvec]
    fvec = [cf1.fvec; cf2.fvec]
    xvec = [cf1.xvec; cf2.xvec]
    initvec = [cf1.initvec; cf2.initvec]
    HyperClosedForm(cf1.func, cf1.arg, evec, rvec, fvec, xvec, initvec)
end

# ------------------------------------------------------------------------------

function (c::CFiniteClosedForm{T})(n::Union{Int, T}) where {T}
    CFiniteClosedForm(c.func, c.arg, c.mvec, c.rvec, [subs(x, c.arg, n) for x in c.xvec], c.initvec, subs(c.instance, c.arg, n))
end

function (c::HyperClosedForm{T})(n::Union{Int, T}) where {T}
    HyperClosedForm(c.func, c.arg, c.evec, c.rvec, c.fvec, [subs(x, c.arg, n) for x in c.xvec], c.initvec, subs(c.instance, c.arg, n))
end

# ------------------------------------------------------------------------------

function reset(c::CFiniteClosedForm{T}) where {T}
    if c.arg in free_symbols(c.instance)
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

function reset(c::HyperClosedForm{T}) where {T}
    if c.arg in free_symbols(c.instance)
        shift = c.instance - c.arg
        factors = [e^shift for e in c.evec]
        xvec = c.xvec .* factors
        evec = c.evec
    else
        xvec = c.xvec .* (c.evec .^ c.instance)
        evec = fill(T(1), length(c.rvec))
    end
    HyperClosedForm(c.func, c.arg, evec, c.rvec, c.fvec, xvec, c.initvec)
end

# ------------------------------------------------------------------------------

init(c::CFiniteClosedForm, d::Dict) = CFiniteClosedForm(c.func, c.arg, c.mvec, c.rvec, [subs(x, d...) for x in c.xvec], c.initvec, c.instance)

init(c::HyperClosedForm, d::Dict) = HyperClosedForm(c.func, c.arg, c.evec, c.rvec, c.fvec, [subs(x, d...) for x in c.xvec], c.initvec, c.instance)

# ------------------------------------------------------------------------------

Base.show(io::IO, cf::Union{CFiniteClosedForm, HyperClosedForm}) = print(io, " $(cf.func)($(cf.instance)) = $(convert(Basic, cf))")

function Base.show(io::IO, ::MIME"text/plain", cf::Union{CFiniteClosedForm, HyperClosedForm})
    summary(io, cf)
    println(io, ":")
    show(io, cf)
end