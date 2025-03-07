using Test

include("../src/smpnp.jl")

@testset "Test SMPNP" begin
    @test multiply(3, 4) == 12
end
