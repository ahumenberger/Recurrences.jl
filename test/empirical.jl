using DataFrames
using MacroTools
using Recurrences

function loop(xs::Vector{Expr}, f::Function; init::Dict, iterations::Int = 10)
    ls = [:(Symbol($(string(v))) => $v) for v in keys(init)]
    is = [:($k = $v) for (k, v) in init]
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

function loop(cs::Vector{<:ClosedForm}, f::Function; init::Dict, iterations::Int = 10)
    xs = [convert(Expr, c) for c in cs]
    xs = [MacroTools.postwalk(x->(x isa Number ? Rational(x) : x), y) for y in xs]
    vs = [Symbol(string(c.func)) for c in cs]
    is = [k in vs ? :($(initvariable(k, 0)) = $v) : :($k = $v) for (k, v) in init]
    ls = [v in vs ? :(Symbol($(string(v))) => $v(it)) : :(Symbol($(string(v))) => $v) for v in keys(init)]
    quote
        let $(is...)
            $(xs...)
            for it in 0:$iterations
                $f(:it=>it, $(ls...))
            end
        end
    end
end

function collecttrace(xs; init, kwargs...)
    df = DataFrame()
    df[:it] = Int[]
    for v in keys(init)
        df[v] = Rational[]
    end
    eval(loop(xs, (x...)->push!(df, Dict(x)); init=init, kwargs...))
    df
end

function cftrace(xs::Vector{Expr}; kwargs...)
    cs = lrs_sequential(xs) |> solve
    collecttrace(cs; kwargs...)
end

function trace(xs::Vector{Expr}; kwargs...)
    vs = Base.unique(Iterators.flatten(map(Recurrences.free_symbols, xs)))
    dict = Dict(v=>Rational(rand(-10:10)) for v in vs)
    ldf = collecttrace(xs; kwargs..., init=dict)
    cdf = cftrace(xs; kwargs..., init=dict)
    ldf, cdf
end