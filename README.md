# Recurrences.jl

A Julia package for solving systems of linear recurrences.

```julia
julia> using Recurrences

julia> lrs = @lrs begin
         x(n+1) = 2*x(n)
         y(n+1) = 1/2*y(n)
       end
2-element LinearRecSystem{SymEngine.Basic} of order 1:
 ⎛ 1  0 ⎞⎛ x(n+1) ⎞ + ⎛ -2     0 ⎞⎛ x(n) ⎞ = ⎛ 0 ⎞
 ⎝ 0  1 ⎠⎝ y(n+1) ⎠   ⎝  0  -1/2 ⎠⎝ y(n) ⎠   ⎝ 0 ⎠

julia> solve(lrs)
2-element Array{CFiniteClosedForm,1}:
  x(n) = 2^n*x00   
  y(n) = 2^(-n)*y00
```
