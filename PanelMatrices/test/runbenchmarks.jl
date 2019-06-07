using PanelMatrices
using BenchmarkTools
using StaticArrays
using StaticNumbers
using LinearAlgebra
using Statistics

# Make sure we're benchmarking BLAS single-threaded
LinearAlgebra.BLAS.set_num_threads(1)

# In case it affects benchmarking
set_zero_subnormals(true)

const y = reshape(collect(1:100), 10, 10)
const x = PanelMatrix(y)
const v = @SMatrix [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16]

println("Benchmarking using Julia", VERSION)

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

for T in (Float32, Float64), M = 7:5:32, N = 7:5:32, K = 7:5:32
    if M*N*K <= 25*25*25
        println("\n$T $M x $N x $K mul!")

        A0 = reshape(collect(T(1):T(M*K)), M, K)
        B0 = reshape(collect(T(1):T(K*N)), K, N)
        C0 = A0*B0

        A = PanelMatrix(A0)
        B = PanelMatrix(B0)
        C = PanelMatrix{T}(undef, (M,N))

        mul!(C, A, B)
        @assert C ≈ C0

        tb = @benchmark mul!($C, $A, $B) samples=10 evals=1000
        println(2*M*N*K/minimum(tb.times), " Gflops")

        println("Julia arrays using $(LinearAlgebra.BLAS.vendor()) (single thread), for comparison")
        tbr = @benchmark mul!($C0, $A0, $B0) samples=100 evals=100
        println(2*M*N*K/minimum(tbr.times), " Gflops")
    end
end
