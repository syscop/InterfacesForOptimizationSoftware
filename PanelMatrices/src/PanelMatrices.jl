module PanelMatrices

using StaticNumbers
using StaticArrays
using UnsafeArrays # maybe we don't need this, but just use the same idea.

#TODO: using UnsafeArrays 

export PanelMatrix

struct ColumnMajorPanel{T,M,N}

"""
A `PanelMatrix` is a panel-major matrix, or a view into a panel-major matrix.
"""
struct PanelMatrix{T,D,M,N,S,O,P,Q,R} <: AbstractMatrix{T}
    data::D # The actual data <: AbstractVector{T}
    size::Tuple{Int, Int} # (rows, columns) in the matrix
    panel_stride::Int # number of panels per column in data

    # The below field are of size zero. Just more convenient than accessing type parameters.
    panel_size::Tuple{StaticInteger{M}, StaticInteger{N}} # (rows, columns) in each panel
    panel_class::Val{S}  # :CM = column major, :RM = row major

    # These may be Int or StaticInteger, depending on tradeoff between compile and runtime.
    pad_first::Tuple{O, P} # (row, column) of first element, within first tile
    pad_last::Tuple{Q, R} # (row, column) of last element, within last tile

    function PanelMatrix(
           data::D
           size::Tuple{Int, Int},
           panel_size::Tupe{StaticInteger{M}, StaticInteger{N}} = default_panelsize(eltype(data))),
           panel_class::<:Val{S} = Val{:CM},
           pad_first::Tuple{O, P} = static.((0, 0)),
           static_pad_last::Val{Z} = Val{false}
           ) where {D<:AbstractArray{T} where T where {T,D,M,N,S,O,P,Q,R}
        @boundscheck begin
            # TODO: checkbounds(data, 1:len)
        end
        pad_last = mod.(panel_size .- offset .- size, panel_size)
        if Z
            pad_last = static.(pad_last)
        end
        (Q, R) = typeof.(pad_last)
        new{T,D,M,N,S,O,P,Q,R}(data, size, panel_size, panel_class, pad_first, pad_last)
   end
end

function PanelMatrix(x::AbstractMatrix, panel_type::Type=default_panel(eltype(x)))
    (rows, columns) = size(x)
    panel_rows=cld(rows, size(panel_type, 1))
    panel_columns=cld(columns, size(panel_type, 2))
    data = zeros(panel_type, (panel_rows, panel_columns)) # TODO: Cache line alignment
    z = PanelMatrix(data, rows, columns)
    #  TODO: implement z .= x
    for i in eachindex(x)
        z[i] = x[i]
    end
    return z
end

default_panelsize(T) = static.((4, 4))

"""
    getpanel(x, i, j)

Get the panel containing element (i,j)
"""
function getpanel(x::PanelMatrix, i, j)
    @boundscheck checkbounds(x, i, j)
    # TODO
end

Base.size(x::PanelMatrix) = x.size
Base.eltype(x::PanelMatrix) = eltype(x.data)

#
# function Base.setindex!(x::PanelMatrix, v, i, j)
#      getpanel(x, i, j)[(i % x.rows_per_panel) + 1 ,(j % x.columns_per_panel) + 1] = v
# end
#
# Base.getindex(x::PanelMatrix, i, j) = getpanel(x, i, j)[(i % x.rows_per_panel) + 1 ,(j % x.columns_per_panel) + 1]

end # module
