using SymPy
using Recurrences: solve

systems = [
(
    @lrs begin
        r(n+1) = r(n) - v(n)
        v(n+1) = v(n) + 2
    end
    ,[
        "-n^2 - n*(v_0 - 1) + r_0",
        "2*n + v_0"
    ]
),
(
    @lrs begin
        r(n+1) = r(n) + u(n)
        u(n+1) = u(n) + 2
    end
    ,[
        "+n^2 + n*(u_0 - 1) + r_0",
        "2*n + u_0"
    ]
)
]

for (lrs, expected) in systems
    exp = simplify.(Sym.(expected))
    expr = expression.(solve(lrs))
    @test Set(expr) == Set(exp)
end
