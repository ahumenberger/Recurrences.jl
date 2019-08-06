# macro lrs(input)
#     T = Basic
#     entries = LinearEntry{T}[]
#     args = T[]
#     @capture(input, begin fields__ end)
#     for ex in fields
#         @capture(ex, lhs_ = rhs_)
#         lhs = sympy.sympify(string(lhs))
#         rhs = sympy.sympify(string(unblock(rhs)))
#         entry, arg = LinearRecEntry(T, lhs - rhs)
#         push!(entries, entry)
#         push!(args, arg)
#     end
#     args = Base.unique(args)
#     @assert length(args) == 1 "More than one function argument, got $(args)"
#     sys = LinearRecSystem(args[1])
#     push!(sys, entries...)
#     return sys
# end

macro lrs(input)
    @capture(input, begin fields__ end)
    exprs = Sym[]
    for ex in fields
        @capture(ex, lhs_ = rhs_)
        lhs = sympy.sympify(string(lhs))
        rhs = sympy.sympify(string(unblock(rhs)))
        push!(exprs, lhs-rhs)
    end
    funcs = Iterators.flatten(function_symbols.(exprs))
    args = Iterators.flatten([f.args for f in funcs]) |> collect
    args = free_symbols(args) |> Base.unique
    @assert length(args) == 1 "More than one function argument, got $(args)"
    lrs(Vector{Expr}(fields), Symbol(string(args[1])))
end

macro lrs_parallel(input)
    @capture(input, begin fields__ end)
    lrs_parallel(Vector{Expr}(fields))
end

macro solve_parallel(input)
    @capture(input, begin fields__ end)
    solve_parallel(Vector{Expr}(fields))
end

macro rec(expr)
    LinearRecEntry(SymPy.Sym(string(expr)))[1]
end

function function_symbols(expr::SymPy.Sym)
    return sympify.(collect(expr.atoms(AppliedUndef)))
end

LinearRecEntry(t::Type{T}, expr::Expr) where {T<:Union{SymPy.Sym, SymEngine.Basic}} = LinearRecEntry(t, sympify(string(expr)))

function LinearRecEntry(::Type{T}, expr::SymPy.Sym) where {T<:Union{SymPy.Sym, SymEngine.Basic}}
    funcs = function_symbols(expr)
    @assert length(funcs) > 0 "Not a recurrence: no functions present"
    args = Iterators.flatten([f.args for f in funcs]) |> collect
    fsyms = free_symbols(args) |> Base.unique
    @assert length(fsyms) == 1 "More (or less) than one variable in the function arguments, got: $(fsyms)"
    
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

function LinearRecEntry(::Type{Var}, ::Type{RAPL}, expr::Expr)
    expr = sympify(string(expr))
    funcs = function_symbols(expr)
    @assert length(funcs) > 0 "Not a recurrence: no functions present"
    args = Iterators.flatten([f.args for f in funcs]) |> collect
    fsyms = free_symbols(args) |> Base.unique
    @assert length(fsyms) == 1 "More (or less) than one variable in the function arguments, got: $(fsyms)"
    
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

    hom = zero(RPoly)
    dicts = [Dict{Var,RPoly}() for _ in 1:maxarg + 1]
    for i in 0:maxarg
        for f in funcs
            fc = f(farg + i)
            co = coeff(expr, f(farg + i))
            if !iszero(co)
                hom += co * fc
                g = mkvar(string(Sym(f)))
                co = mkpoly(string(co))
                dicts[i+1][g] = co
            end
        end
    end

    @info "" string(-(expr - hom))
    inhom = mkpoly(string(-(expr - hom)))
    farg = mkvar(string(farg))

    (coeffs = dicts, inhom = inhom), farg
end

free_symbols(x::Union{SymPy.Sym, Vector{SymPy.Sym}}) = SymPy.free_symbols(x)
free_symbols(x::Union{SymEngine.Basic, Vector{SymEngine.Basic}}) = SymEngine.free_symbols(x)

coeff(x::SymPy.Sym, b::SymPy.Sym) = x.coeff(b)

subs(x::SymPy.Sym, args...) = SymPy.subs(x, args...)
subs(x::SymEngine.Basic, args...) = SymEngine.subs(x, args...)

subs(ex::AbstractArray, args...; kwargs...) = map(u -> subs(u, args...; kwargs...), ex)
# subs(ex::AbstractArray{SymEngine.Basic}, args...; kwargs...) = map(u -> subs(u, args...; kwargs...), ex)

Base.inv(m::Matrix{SymEngine.Basic}) = convert(Matrix{SymEngine.Basic}, inv(convert(Matrix{Sym}, m)))

expand(x::Sym) = SymPy.expand(x)
expand(x::Basic) = SymEngine.expand(x)

solve(x::Vector{Sym}, y::Vector{Sym}) = SymPy.solve(x, y)
solve(x::Vector{Basic}, y::Vector{Basic}) = convert(Dict{Basic,Basic}, solve(convert.(Sym, x), convert.(Sym, y)))

coeff(x::Basic, y::Basic) = SymEngine.coeff(x, y)

function coeffs(p::Sym, x::Sym)
    if iszero(p)
        return []
    end
    SymPy.coeffs(p, x)
end
function coeffs(p::Basic, x::Basic)
    convert.(Basic, coeffs(convert(Sym, p), convert(Sym, x)))
end

Base.:\(A::Matrix{Basic}, b::Vector{Basic}) = convert(Vector{Basic}, convert(Matrix{Sym}, A) \ convert(Vector{Sym}, b))