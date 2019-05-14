using PanelMatrices
using Test
using StaticArrays
using StaticNumbers

@testset "PanelMatrices.jl" begin
    y = reshape(collect(1:100), 10, 10)
    x = PanelMatrix(y)
    @test x[1,1] == y[1,1]
    @test all(x .== y)

    v = PanelMatrices.panel_view(x, 1, 1)
    @test all(v .== [1 11 21 31; 2 12 22 32; 3 13 23 33; 4 14 24 34])
    @test all(PanelMatrices.panel_view(x, 3, 1) .== [9 19 29 39; 10 20 30 40])
    @test all(PanelMatrices.get_panel(x, 3, 1) .== [9 19 29 39; 10 20 30 40])
    @test all(PanelMatrices.panel_view(x, 1, 3) .== [81 91; 82 92; 83 93; 84 94])
    @test all(PanelMatrices.get_panel(x, 1, 3) .== [81 91; 82 92; 83 93; 84 94])
    @test all(PanelMatrices.panel_view(x, 3, 3) .== [89 99; 90 100])
    @test all(PanelMatrices.get_panel(x, 3, 3) .== [89 99; 90 100])
    @test size(PanelMatrices.full_panel_view(x, 3, 3)) == (4, 4)

    p = PanelMatrices.get_full_panel(x, 1, 1)
    @test p === @SMatrix [1 11 21 31; 2 12 22 32; 3 13 23 33; 4 14 24 34]

    z = reshape(101:116, 4, 4)
    PanelMatrices.set_full_panel!(x, z, 2, 2)
    @test x[5,5] == z[1,1]
    @test all([x[i,j] for i ∈ 5:8, j ∈ 5:8] .== z)

    # TODO: Write tests that ensure that there is no allocation in critical routines.
end
