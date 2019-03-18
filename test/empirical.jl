using DataFrames
using MacroTools
using Recurrences

function free_symbols(ex::Expr)
    ls = Symbol[]
    MacroTools.postwalk(x -> x isa Symbol && Base.isidentifier(x) ? push!(ls, x) : x, ex)
    Base.unique(ls)
end

getval(x::Module, s::Symbol) = getfield(x, s)
getval(x::Dict, s::Symbol) = x[s]

function loop(xs::Vector{Expr}, f::Function; iterations::Int = 10, initscope = nothing)
    vs = Base.unique(Iterators.flatten(map(free_symbols, xs)))
    ls = [:(Symbol($(string(v))) => $v) for v in vs]
    if initscope == nothing
        is = [:($v = $(Rational(rand(-10:10)))) for v in vs]
    else
        is = [:($v = $(getval(initscope, v))) for v in vs]
    end
    quote
        let $(is...)
            $f(:it=>0, $(ls...))
            for it in 1:$iterations
                $(xs...)
                $f(:it=>it, $(ls...))
            end
        end
    end
end

function loop(cs::Vector{<:ClosedForm}, f::Function; iterations::Int = 10, initscope = nothing)
    xs = [convert(Expr, c) for c in cs]
    vs = [Symbol(string(c.func)) for c in cs]
    ls = [:(Symbol($(string(v))) => $v(it)) for v in vs]
    # ls0 = [:(Symbol($(string(v))) => $(initvariable(v, 0))) for v in vs]
    if initscope == nothing
        is = [:($(initvariable(v, 0)) = $(Rational(rand(-10:10)))) for v in vs]
    else
        is = [:($(initvariable(v, 0)) = $(getval(initscope, v))) for v in vs]
    end
    quote
        let $(is...)
            $(xs...)
            for it in 0:$iterations
                $f(:it=>it, $(ls...))
            end
        end
    end
end

function looptrace(xs; kwargs...)
    df = DataFrame()
    if xs isa Vector{Expr}
        vs = Base.unique(Iterators.flatten(map(free_symbols, xs)))
    else
        vs = [Symbol(string(x.func)) for x in xs]
    end
    df[:it] = Int[]
    for v in vs
        df[v] = Rational[]
    end
    eval(loop(xs, (x...)->push!(df, Dict(x)); kwargs...))
    df
end

function cftrace(xs::Vector{Expr}; kwargs...)
    cs = lrs_sequential(xs) |> solve
    looptrace(cs; kwargs...)
end

function trace(xs::Vector{Expr}; kwargs...)
    vs = Base.unique(Iterators.flatten(map(free_symbols, xs)))
    dict = Dict(v=>Rational(rand(-10:10)) for v in vs)
    @info "" dict
    ldf = looptrace(xs; kwargs..., initscope=dict)
    cdf = cftrace(xs; kwargs..., initscope=dict)
    ldf, cdf
end