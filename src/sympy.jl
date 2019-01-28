export @rec, @lrs

macro lrs(input)
    entries = LinearEntry{SymPy.Sym}[]
    args = SymPy.Sym[]
    @capture(input, begin fields__ end)
    for ex in fields
        @capture(ex, lhs_ = rhs_)
        lhs = SymPy.Sym(string(lhs))
        rhs = SymPy.Sym(string(unblock(rhs)))
        entry, arg = LinearRecEntry(lhs - rhs)
        push!(entries, entry)
        push!(args, arg)
    end
    args = Base.unique(args)
    @assert length(args) == 1 "More than one function argument, got $(args)"
    sys = LinearRecSystem(args[1])
    push!(sys, entries...)
    sys
end

macro rec(expr)
    LinearRecEntry(SymPy.Sym(string(expr)))[1]
end

function function_symbols(expr::SymPy.Sym)
    return Sym.(collect(atoms(expr, AppliedUndef)))
end

function LinearRecEntry(expr::SymPy.Sym)
    funcs = function_symbols(expr)
    @assert length(funcs) > 0 "Not a recurrence: no functions present"
    args = Iterators.flatten(SymPy.args.(funcs)) |> collect
    fsyms = free_symbols(args) |> Base.unique
    @assert length(fsyms) == 1 "More (or less) than one variable in the function arguments, got: $(fsyms)"
    
    farg = fsyms[1]
    # remove functions which are not of the form x(n+1)
    funcs = func.(filter(x -> has(x, farg), funcs)) |> Base.unique
    args = filter(x -> has(x, farg), args)
    @assert !(string(farg) in string.(funcs)) "Ambiguous symbol: $(farg) is a function and an argument"

    args = args .- farg
    minarg = convert(Int, minimum(args))
    if !iszero(minarg)
        shifted = subs(expr, farg => farg - minarg)
        return LinearRecEntry(shifted)
    end
    maxarg = convert(Int, maximum(args))

    hom = SymPy.Sym(0)
    dicts = [Dict{SymPy.Sym, SymPy.Sym}() for _ in 1:maxarg + 1]
    for i in 0:maxarg
        for f in funcs
            fc = f(farg + i)
            co = coeff(expr, f(farg + i))
            if !iszero(co)
                hom += co * fc
                dicts[i+1][f] = co
            end
        end
    end

    (coeffs = dicts, inhom = -(expr - hom)), farg
end