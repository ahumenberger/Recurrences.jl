using BenchmarkTools
using Recurrences

include("../test/examples.jl")

const bench = BenchmarkGroup()
bench["cfinite"] = BenchmarkGroup()

for bsym in branches
    b = eval(bsym)
    bench["cfinite"][string(bsym)] = @benchmarkable lrs_sequential($b)
end

tune!(bench)
results = run(bench)

showall(results)