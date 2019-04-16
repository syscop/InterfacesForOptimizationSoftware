using PanelMatrices
using BenchmarkTools
using StaticArrays

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
