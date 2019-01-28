export Recurrence, LinearRecurrence, CFiniteRecurrence, CFiniteClosedForm

abstract type Recurrence end
abstract type LinearRecurrence <: Recurrence end

struct CFiniteRecurrence{T} <: LinearRecurrence
    func::T
    arg::T
    coeffs::Vector{T}
    inhom::T

    function CFiniteRecurrence{T}(func::T, arg::T, coeffs::Vector{T}, inhom::T = T(0)) where {T}
        if !all(is_constant.(coeffs))
            error("Not a C-finite recurrence.")
        end
        new{T}(func, arg, coeffs, inhom)
    end
end

coeffs(r::CFiniteRecurrence) = r.coeffs
order(r::CFiniteRecurrence) = length(r.coeffs) - 1

struct CFiniteClosedForm{T}
    func::T
    arg::T
    mvec::Vector{T}
    rvec::Vector{T}
    xvec::Vector{T}
    initvec::Vector{T}
end

function mroots(poly::Polynomials.Poly{T}) where {T}
    roots = Polynomials.roots(poly)
    Dict([(r, count(x -> x==r, roots)) for r in Base.unique(roots)])
end

function mroots(p::Polynomials.Poly{SymPy.Sym})
    T = SymPy.Sym
    length(p) == 0 && return zeros(T, 0)

    num_leading_zeros = 0
    while p[num_leading_zeros] ≈ zero(T)
        if num_leading_zeros == length(p)-1
            return zeros(R, 0)
        end
        num_leading_zeros += 1
    end
    num_trailing_zeros = 0
    while p[end - num_trailing_zeros] ≈ zero(T)
        num_trailing_zeros += 1
    end
    n = lastindex(p)-(num_leading_zeros + num_trailing_zeros)
    n < 1 && return zeros(T, length(p) - num_trailing_zeros - 1)

    companion = diagm(-1 => ones(T, n-1))
    an = p[end-num_trailing_zeros]
    companion[1,:] = -p[(end-num_trailing_zeros-1):-1:num_leading_zeros] / an

    companion[:eigenvals]()
end


function closedform(rec::CFiniteRecurrence{T}) where {T}

    # TODO: allow inhomogeneous recurrences?
    roots = Polynomials.Poly(coeffs(rec)) |> mroots
    # @info "Roots" roots

    size = order(rec)
    mvec = [T(i) for (_, m) in roots for i in 0:m - 1] # multiplicities
    rvec = [z for (z, m) in roots for _ in 0:m - 1] # roots

    A = [i^m * r^i for i in 0:size-1, (r, m) in zip(rvec, mvec)]
    b = [T(string(string(rec.func), i)) for i in 0:size - 1] 
    # @info "Ansatz" A b A\b
    CFiniteClosedForm(rec.func, rec.arg, mvec, rvec, A \ b, b)
end

function Base.show(io::IO, cf::CFiniteClosedForm)
    vec = [cf.arg^m * r^cf.arg for (r, m) in zip(cf.rvec, cf.mvec)]
    res = simplify(transpose(vec) * cf.xvec)
    println(io, "$(typeof(cf)):")
    print(io, " $(cf.func)($(cf.arg)) = $(string(res))")
end