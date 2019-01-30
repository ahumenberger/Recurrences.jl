using SymPy
using Recurrences: solve

lrs_fermat_b1 = @lrs begin
    r(n+1) = r(n) - v(n)
    v(n+1) = v(n) + 2
end

lrs_fermat_b1_res = lrs_fermat_b1 |> solve

@info lrs_fermat_b1

# lrs_fermat_b1_exp = [
#     Sym("-n^2 - n*(v_0 - 1) + r_0")
#     Sym("2*n + v_0")
# ]

lrs_fermat_b2 = @lrs begin
    r(n+1) = r(n) + v(n)
    u(n+1) = u(n) + 2
end



@test expression.(solve(lrs_fermat_b1)) == nothing