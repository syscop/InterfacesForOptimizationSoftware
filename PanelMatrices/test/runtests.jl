using PanelMatrices
using Test
using StaticArrays
using StaticNumbers

@testset "PanelMatrices.jl" begin
    y = reshape(1:100, 10, 10)
    x = PanelMatrix(y)
    println(x)
    println(x.data)
    println(y)
    @test x[1,1] == y[1,1]
    @test all(x .== y)

    # TODO: Write tests that ensure that there is no allocation
end
