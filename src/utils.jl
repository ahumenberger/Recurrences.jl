
function variables(::Type{T}; n::Int = 1, unique::Bool = true) where {T}
    if unique
        global var_count += n
        varcnt = var_count
    else
        varcnt = n
    end
    if n == 1
        return var(T, "ω$(varcnt)")
    end
    return [var(T, "ω$i") for i in varcnt-n+1:varcnt]
end

var(::Type{SymPy.Sym}, s::String) = SymPy.Sym(s)
var(::Type{SymEngine.Basic}, s::String) = SymEngine.symbols(s)

function summands(expr::Sym)
    expr = expand(expr)
    @debug "" expr typeof(expr)
    if SymPy.funcname(expr) == "Add"
        return args(expr)
    end
    [expr]
end
summands(x::SymEngine.Basic) = convert.(Basic, summands(convert(Sym, x)))

lcm2(n::Sym, rest::Sym...) = SymPy.lcm(n, lcm2(rest...))
lcm2(x::Basic...) = convert(Basic, lcm2(convert.(Sym, x)...))
lcm2() = 1

function pascal(n; alt = false)
    f = alt ? -1 : 1
    entries = [1]
    for k in 0:n-1
        append!(entries, f * entries[k+1] * (n-k) / (k+1))
    end
    entries
end

function free_symbols(ex::Expr)
    ls = Symbol[]
    MacroTools.postwalk(x -> x isa Symbol && Base.isidentifier(x) ? push!(ls, x) : x, ex)
    Base.unique(ls)
end

initvar(v::T, i::Union{T, Int64}=0) where {T} = T("$(string(v))$(i)$(i)")