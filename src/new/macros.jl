replace_post(ex, s, s′) = postwalk(x -> x == s ? s′ : x, ex)

function LinearRecEntry(expr::Expr)
    CoeffT = QQ
    _funcs = Symbol[]
    _args = RExpr[]
    function_walk(expr) do f, x
        @assert length(x) == 1
        push!(_funcs, f)
        push!(_args, x[1])
    end
    _funcs = Base.unique(_funcs)
    @assert length(_funcs) > 0 "Not a recurrence: no functions present"
    argsyms = Base.unique(Iterators.flatten(symbols(a) for a in _args))
    @assert length(argsyms) == 1 "More (or less) than one variable in the function arguments, got: $(argsyms)"

    _argsym = argsyms[1]
    @assert !(_argsym in _funcs) "Ambiguous symbol: $(_argsym) is a function and an argument"

    R, argsym = PolynomialRing(CoeffT, string(_argsym))

    s_expr = SymEngine.Basic(expr)
    s_argsym = SymEngine.Basic(_argsym)
    s_args = map(x->(SymEngine.Basic(x) - s_argsym), _args)
    # minarg = minimum(convert(Rational, x) for x in s_args)
    minarg = convert(Int, minimum(s_args))
    if !iszero(minarg)
        _expr = replace_post(expr, _argsym, :($_argsym - $minarg))
        return LinearRecEntry(_expr)
    end
    maxarg = convert(Int, maximum(s_args))

    hom = zero(SymEngine.Basic)
    dicts = [Dict{Symbol,FracElem}() for _ in 1:maxarg + 1]
    for i in 0:maxarg
        for f in _funcs
            fc = SymEngine.SymFunction(f)(s_argsym + i)
            co = SymEngine.coeff(s_expr, fc)
            if !iszero(co)
                hom += co * fc
                # call R twice to deal with integers
                num = numerator(co)
                den = denominator(co)
                cnum = R(R(convert(Expr, num)))
                cden = R(R(convert(Expr, den)))
                dicts[i+1][f] = cnum // cden
            end
        end
    end

    inhom = R(R(convert(Expr, SymEngine.expand(-(s_expr - hom)))))//R(1)

    (coeffs = dicts, inhom = inhom), argsym
end

macro lrs(input)
    @capture(input, begin fields__ end)
    lrs(Vector{Expr}(fields))
end

function entries(ex::Expr...)
    _entries = []
    _args = []
    _syms = Symbol[]
    for x in ex
        if @capture(x, l_ = r_)
            x = :($(unblock(l)) - $(unblock(r)))
            @assert @capture(l, v_(a_))
            push!(_syms, v)
        end
        entry, arg = LinearRecEntry(x)
        push!(_entries, entry)
        push!(_args, arg)
    end
    _args = Base.unique(_args)
    @assert length(_args) == 1
    _syms, _entries, _args[1]
end

function lrs(exprs::Vector{Expr})
    syms, ls, arg = entries(exprs...)
    lrs = LinearRecSystem(arg, syms)
    push!(lrs, ls...)
    lrs
end