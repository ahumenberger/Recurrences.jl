var_count = 0

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

var(::Type{Symbol}, s::String) = Symbol(s)

function pascal(n; alt = false)
    f = alt ? -1 : 1
    entries = [1]
    for k in 0:n-1
        append!(entries, f * entries[k+1] * (n-k) / (k+1))
    end
    entries
end

initvar(v::T, i::Union{T, Int64}=0) where {T} = T("$(string(v))$(i)$(i)")