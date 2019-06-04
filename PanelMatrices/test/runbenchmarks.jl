using PanelMatrices
using BenchmarkTools
using StaticArrays
using StaticNumbers
using LinearAlgebra

const y = reshape(collect(1:100), 10, 10)
const x = PanelMatrix(y)
const v = @SMatrix [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]

println("Benchmarking using Julia ", VERSION)

let
    println("Setting each individual element")
    function f(x, s)
        for i ∈ eachindex(x)
            @inbounds x[i] = s
        end
    end
    @btime $f($x, 0)
end

let
    println("Setting each panel")
    function f(x, v)
        for i ∈ 1:3, j ∈ 1:3
            @inbounds PanelMatrices.set_full_panel!(x, v, i, j)
        end
    end
    @btime $f($x, $v)
end

let
    println("16x16x16 matmul")
    A0 = reshape(collect(1.0:256.0), 16, 16)
    B0 = reshape(collect(1.0:256.0), 16, 16)
    C0 = A0*B0

    A = PanelMatrix(A0, static.((4,4)))
    B = PanelMatrix(B0, static.((4,4)))
    C = PanelMatrix{Float64}(undef, (16,16), static.((4,4)))

    @btime mul!($C, $A, $B)

    # println("16x16x16 matmul (Julia arrays, for comparison)")
    # @btime mul!($C0, $A0, $B0)
    #
    # println("16x16x16 matmul (MMatrix, for comparison)")
    # As = MMatrix{size(A0)...}(A0)
    # Bs = MMatrix{size(B0)...}(B0)
    # Cs = MMatrix{size(B0)...}(C0)
    # @btime mul!($Cs, $As, $Bs)
end
