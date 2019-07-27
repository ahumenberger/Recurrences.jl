replace_post(ex, s, s′) = MacroTools.postwalk(x -> x == s ? s′ : x, ex)
gensym_unhashed(s::Symbol) = Symbol(replace(string(gensym(s)), "#"=>""))

const RExpr = Union{Expr,Symbol,Number}

function split_assign(xs::Vector{Expr})
    ls = Symbol[]
    rs = RExpr[]
    for x in xs
        @capture(x, l_ = r_)
        push!(ls, unblock(l))
        push!(rs, unblock(r))
    end
    ls, rs
end

function pushexpr!(lrs::LinearRecSystem, ls::Vector{Expr}, rs::Vector{RExpr})
    for (l, r) in zip(ls, rs)
        expr = :($l - $r)
        entry, _ = LinearRecEntry(Basic, expr)
        push!(lrs, entry)
    end
    lrs
end

function pushexpr!(lrs::LinearRecSystem, ex::Expr...)
    for x in ex
        if @capture(x, l_ = r_)
            x = :($l - $r)
        end
        entry, _ = LinearRecEntry(Basic, x)
        push!(lrs, entry)
    end
    lrs
end

function lrs_sequential(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lhss, rhss = split_assign(exprs)
    lrs = LinearRecSystem(Basic(lc), map(Basic, lhss))
    for (i,v) in enumerate(lhss)
        rhss = RExpr[replace_post(x, v, (i < j ? :($v($lc+1)) : :($v($lc)))) for (j,x) in enumerate(rhss)]
    end
    lhss = [Expr(:call, v, Expr(:call, :+, lc, 1)) for v in lhss]
    pushexpr!(lrs, lhss, rhss)
end

# ------------------------------------------------------------------------------

function solve_sequential(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lhss, rhss = split_assign(exprs)
    for (i,v) in enumerate(lhss)
        rhss = RExpr[replace_post(x, v, (i < j ? :($v($lc+1)) : :($v($lc)))) for (j,x) in enumerate(rhss)]
    end
    sys = [:($l($lc+1) - $r) for (l, r) in zip(lhss, rhss)]
    @debug "Solve sequential" sys
    solve(sys, lhss, lc)
end

function solve(exprs::Vector{Expr}, vars::Vector{Symbol}, lc::Symbol = gensym_unhashed(:n))
    lrs = LinearRecSystem(Basic(lc))
    nonlinear = Expr[]
    for expr in exprs
        entry, _ = LinearRecEntry(Basic, expr)
        if entry == nothing
            push!(nonlinear, expr)
        else
            push!(lrs, entry)
        end
    end
    if !isempty(lrs)
        @debug "Nonempty LRS" lrs
        # for lin in linear

        # end
        # push!(lrs, linear...)
        cfs = solve(lrs)
        if isempty(nonlinear)
            return cfs
        end
        newsys = Expr[]
        newvars = Symbol[]
        for nl in nonlinear
            fsyms = function_symbols(nl)
            for fsym in fsyms
                f = Basic(fsym.args[1])
                a = Basic(fsym.args[2])
                idx = findfirst(x->x.func == f, cfs)
                if idx != nothing
                    repl = rhs(Expr, cfs[idx](Basic(a)))
                    nl = replace_post(nl, fsym, repl)
                else
                    push!(newvars, fsym.head)
                end
            end
            push!(newsys, nl)
        end
        newcfs =  solve(newsys, newvars, lc)
        return [cfs; newcfs]
    end
    @error("Cannot solve recurrence system.")
end

# ------------------------------------------------------------------------------

function lrs_parallel(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lhss, rhss = split_assign(exprs)
    lrs = LinearRecSystem(Basic(lc), map(Basic, lhss))
    for (i, rhs) in enumerate(rhss)
        rhss[i] = MacroTools.postwalk(x -> x isa Symbol && x in lhss ? :($x($lc)) : x, rhs)
    end
    lhss = [Expr(:call, v, Expr(:call, :+, lc, 1)) for v in lhss]
    pushexpr!(lrs, lhss, rhss)
end

function lrs(exprs::Vector{Expr}, lc::Symbol = gensym_unhashed(:n))
    lrs = LinearRecSystem(Basic(lc))
    pushexpr!(lrs, exprs...)
end

macro lrs(input)
    T = Basic
    entries = LinearEntry{T}[]
    args = T[]
    @capture(input, begin fields__ end)
    for ex in fields
        @capture(ex, lhs_ = rhs_)
        lhs = sympify(string(lhs))
        rhs = sympify(string(unblock(rhs)))
        entry, arg = LinearRecEntry(T, lhs - rhs)
        push!(entries, entry)
        push!(args, arg)
    end
    args = Base.unique(args)
    @assert length(args) == 1 "More than one function argument, got $(args)"
    sys = LinearRecSystem(args[1])
    push!(sys, entries...)
    return sys
end

macro rec(expr)
    LinearRecEntry(SymPy.Sym(string(expr)))[1]
end

function function_symbols(expr::SymPy.Sym)
    return sympify.(collect(expr.atoms(AppliedUndef)))
end

deglist(x, vars) = convert(Tuple, degree_list(x, vars...).x)

islinear(x, vars) = isempty(vars) || all(sum(deglist(a, vars)) <= 1 for a in args(expand(x)))

LinearRecEntry(t::Type{T}, expr::Expr) where {T<:Union{SymPy.Sym, SymEngine.Basic}} = LinearRecEntry(t, Sym(string(expr)))

function LinearRecEntry(::Type{T}, expr::SymPy.Sym) where {T<:Union{SymPy.Sym, SymEngine.Basic}}
    funcs = function_symbols(expr)
    @assert length(funcs) > 0 "Not a recurrence: no functions present"
    args = Iterators.flatten([f.args for f in funcs]) |> collect
    fsyms = free_symbols(args) |> Base.unique
    @assert length(fsyms) == 1 "More (or less) than one variable in the function arguments, got: $(fsyms)"
    # if !islinear(expr, funcs)
    if !sympy.Poly(expr, funcs...).is_linear
        return nothing, nothing
    end

    farg = fsyms[1]
    # remove functions which are not of the form x(n+1)
    funcs = [x.func for x in filter(x -> x.has(farg), funcs)] |> Base.unique
    args = filter(x -> x.has(farg), args)
    @assert !(string(farg) in string.(funcs)) "Ambiguous symbol: $(farg) is a function and an argument"

    args = args .- farg
    minarg = convert(Int, minimum(args))
    # if !iszero(minarg)
    #     shifted = subs(expr, farg => farg - minarg)
    #     return LinearRecEntry(T, shifted)
    # end
    maxarg = convert(Int, maximum(args))

    hom = zero(SymPy.Sym)
    dicts = [Dict{T,T}() for _ in 1:maxarg + 1]
    for i in 0:maxarg
        for f in funcs
            fc = f(farg + i)
            co = coeff(expr, f(farg + i))
            if !iszero(co)
                hom += co * fc
                g = var(T, string(Sym(f)))
                if T == SymEngine.Basic
                    co = convert(SymEngine.Basic, string(co))
                end
                dicts[i+1][g] = co
            end
        end
    end

    inhom = -(expr - hom)
    if T == SymEngine.Basic
        inhom = convert(SymEngine.Basic, string(inhom))
        farg = convert(SymEngine.Basic, string(farg))
    end


    (coeffs = dicts, inhom = inhom), farg
end