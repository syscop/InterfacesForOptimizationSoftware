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

for T in (Float32, Float64), n in (16, 32, 64)
    println("\n$T $n x $n x $n mul!")
    A0 = reshape(collect(one(T):T(n*n)), n, n)
    B0 = reshape(collect(one(T):T(n*n)), n, n)
    C0 = A0*B0

    A = PanelMatrix(A0)
    B = PanelMatrix(B0)
    C = PanelMatrix{T}(undef, (n,n))

    @btime mul!($C, $A, $B)

    println("Julia arrays using $(LinearAlgebra.BLAS.vendor()), for comparison")
    @btime mul!($C0, $A0, $B0)

    # println("StaticArrays.MMatrix, for comparison")
    # As = MMatrix{size(A0)...}(A0)
    # Bs = MMatrix{size(B0)...}(B0)
    # Cs = MMatrix{size(B0)...}(C0)
    # @btime mul!($Cs, $As, $Bs)
end
