module PanelMatrices

using StaticNumbers

"""
PanelMatrix is panel-major, with column-major ordering of panels and within panels.
"""
struct PanelMatrix{T,V,M,N,Q,R,O,P,U,V} <: AbstractMatrix{T}
    data::V # (view on) the actual memory
    eltype::Type{T} # type of the elements
    rows::StaticInteger{M} # rows in this matrix
    cols::StaticInteger{N} # columns in this matrix
    rows_per_panel::StaticInteger{Q} # rows in a panel (rows in memory)
    cols_per_panel::StaticInteger{R} # columns in a panel (columns in memory)
    offset_row::StaticInteger{O} # starting row of the matrix (in memory)
    offset_col::StaticInteger{P} # starting column of the matrix (in memory)
    panel_rows::StaticInteger{U} # number of rows of panels
    panel_cols::StaticInteger{V} # number of columns of panels
    function PanelMatrix(data::AbstractArray{T}, rows, cols, rows_per_panel, cols_per_panel,
                         offset_row=0, offset_col=0,
                         panel_rows=ceil(Int, (rows+offset_row)//rows_per_panel),
                         panel_cols=ceil(Int, (cols+offset_col)//cols_per_panel),
                         ) where T
        @boundscheck length(data) >= panel_rows*rows_per_panel*panel_cols*cols_per_panel || error("Too little data!")
        new{data, T, }() #TODO
    end
end

"""
Panel has column-major ordering of elements.
"""
struct Panel{T,V,M,N,Q,R,O,P} <: AbstractMatrix{T}
    data::V # view on the actual memory
    eltype::Type{T} # type of the elements
    rows::StaticInteger{M} # rows of the view
    cols::StaticInteger{N} # columns of the view
    rows_per_panel::StaticInteger{Q} # rows in this panel (rows in memory)
    cols_per_panel::StaticInteger{R} # columns in this panel (columns in memory)
    offset_row::StaticInteger{O} # starting row of the view
    offset_col::StaticInteger{P} # starting column of the view
    function Panel(data::AbstractArray{T}, rows, cols, rows_per_panel, cols_per_panel,
                         offset_row=0, offset_col=0) where T
        @boundscheck length(data) >= rows_per_panel*cols_per_panel || error("Too little data!")
        new{data, T, }() #TODO
    end
end


# TODO: We need paneling!
getindex(x::Panel, i, j) = x.data[(j+x.offset_col)*x.cols_+(i+x.offset_row)*x.rows_stride]

end # module
