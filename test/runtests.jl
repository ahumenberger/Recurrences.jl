using Test
using Recurrences

include("examples.jl")
include("empirical.jl")

for x in keys(cf)
    @info "$(string(x))"
    a, b = trace(eval(x))
    @info "" a b
    @test a == b
    # @test ==(trace(eval(x))...)
end

# for x in keys(cf)
#     lrs = lrs_sequential(eval(x), :n)
#     @test solve(lrs) == cf[x]
# end