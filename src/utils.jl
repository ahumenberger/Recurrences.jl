
function variables(::Type{Sym}; n::Int = 1, unique::Bool = true)
    if unique
        global var_count += n
        varcnt = var_count
    else
        varcnt = n
    end
    if n == 1
        return Sym("ω$(varcnt)")
    end
    return [Sym("ω$i") for i in varcnt-n+1:varcnt]
end

function summands(expr::Sym)
    expr = expand(expr)
    @debug "" expr typeof(expr)
    if SymPy.funcname(expr) == "Add"
        return args(expr)
    end
    [expr]
end

lcm2(n::Sym, rest::Sym...) = SymPy.lcm(n, lcm2(rest...))
lcm2() = 1

function pascal(n; alt = false)
    f = alt ? -1 : 1
    entries = [1]
    for k in 0:n-1
        append!(entries, f * entries[k+1] * (n-k) / (k+1))
    end
    entries
end