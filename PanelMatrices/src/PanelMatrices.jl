module PanelMatrices

using StaticNumbers

struct PanelMatrix{T,V,M,N,O,P,Q,R,S,T}
    eltype::Type{T}
    data::V
    rows::StaticInteger{M}
    cols::StaticInteger{N}
    offset_row::StaticInteger{O}
    offset_col::StaticInteger{P}
    panel_rows::StaticInteger{Q}
    panel_cols::StaticInteger{R}
    rows_stride::StaticInteger{S}
    cols_stride::StaticInteger{T}
end

end # module
