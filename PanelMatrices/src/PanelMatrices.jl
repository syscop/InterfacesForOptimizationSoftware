module PanelMatrices

using StaticNumbers
using StaticArrays

export PanelMatrix

"""
A `PanelMatrix` is a panel-major matrix, or a view into a panel-major matrix.

It contains a `.data` field which holds the data

Tiles can be, for example, of type `StaticMatrix`.
"""
struct PanelMatrix{T,D,M,N,U,V,O,P,Q,R} <: AbstractMatrix{T}
    data::D # the actual data (must be an AbstractMatrix{S} of panels)
    element_type::Type{T} # type of the elements
    rows::StaticInteger{M} # rows in this matrix
    cols::StaticInteger{N} # columns in this matrix
    rows_per_panel::StaticInteger{U} # number of rows in a panel
    cols_per_panel::StaticInteger{V} # number of columns in a panel
    offset_row::StaticInteger{O} # starting row of the matrix
    offset_col::StaticInteger{P} # starting column of the matrix
    panel_rows::StaticInteger{Q} # number of rows of panels (can be larger than necessary)
    panel_cols::StaticInteger{R} # number of columns of panels (can be larger than necessary)
end

"""
A `Panel` is column major
"""
struct Panel{T,D,M,N,S,U,V,O,P} <: AbstractMatrix{T}
    data::D # the actual data (must be an AbstractMatrix{S} of panels)
    element_type::Type{T} # type of the elements
    rows::StaticInteger{M} # rows in this matrix
    cols::StaticInteger{N} # columns in this matrix
    panel_type::Type{S} # Type of each panel
    rows_per_panel::StaticInteger{U} # number of rows in a panel
    cols_per_panel::StaticInteger{V} # number of columns in a panel
    offset_row::StaticInteger{O} # starting row of the matrix
    offset_col::StaticInteger{P} # starting column of the matrix
end

function PanelMatrix(data::AbstractArray{T}, rows, cols,
                     offset_row=0, offset_col=0,
                     panel_rows=cld(rows+offset_row, size(S,1)),
                     panel_cols=cld(cols+offset_col, size(S,2)),
                     ) where S
    (rows_per_panel, cols_per_panel) = size(S)
    @boundscheck begin
        length(data) >= panel_rows*panel_cols*rows_per_panel*cols_per_panel || error("Too little data!") # TODO Create proper exceptopns as constants
        offset_row + rows <= panel_rows*rows_per_panel || error("Too fews rows of panels")
        offset_col + cols <= panel_cols*cols_per_panel || error("Too fews columns of panels")
        all((rows, cols, offset_row, offset_col) .>= 0) || error("Negative sizes/offsets are not allowed.")
    end
    return PanelMatrix(data, eltype(S), static(Int(rows)), static(Int(cols)), S,
        static(Int(rows_per_panel)), static(Int(cols_per_panel)),
        static(Int(offset_row)), static(Int(offset_col)),
        static(Int(panel_rows)), static(Int(panel_cols)))
end

function PanelMatrix(x::AbstractMatrix, panel_type::Type=default_panel(eltype(x)))
    (rows, cols) = size(x)
    panel_rows=cld(rows, size(panel_type, 1))
    panel_cols=cld(cols, size(panel_type, 2))
    data = zeros(panel_type, (panel_rows, panel_cols)) # TODO: Cache line alignment
    z = PanelMatrix(data, rows, cols)
    #  TODO: implement z .= x
    for i in eachindex(x)
        z[i] = x[i]
    end
    return z
end

default_panel(T) = MMatrix{4, 4, T, 16}

"""
    getpanel(x, i, j)

Get the panel containing element (i,j)
"""
function getpanel(x::PanelMatrix, i, j)
    @boundscheck checkbounds(x, i, j)
    # TODO inbounds
    #x.data[fld(i+x.offset_row, x.rows_per_panel)+1, fld(j+x.offset_col, x.cols_per_panel)+1]
    Panel(view(x.data, ))
end

#TODO: Benchmark whether fixed-length views allocate.

Base.size(x::PanelMatrix) = (x.rows, x.cols)

Base.eltype(x::PanelMatrix) = x.element_type

function Base.setindex!(x::PanelMatrix, v, i, j)
     getpanel(x, i, j)[(i % x.rows_per_panel) + 1 ,(j % x.cols_per_panel) + 1] = v
end

Base.getindex(x::PanelMatrix, i, j) = getpanel(x, i, j)[(i % x.rows_per_panel) + 1 ,(j % x.cols_per_panel) + 1]

end # module
