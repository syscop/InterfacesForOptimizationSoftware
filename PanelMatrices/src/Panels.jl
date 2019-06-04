"""
A `Panel` represents one of the panels in a Panel matrix.

It is used like a view into the matrix.
"""
struct Panel{T,D,P1,P2,C,F1,F2,L1,L2} <: AbstractMatrix{T}
    data::D # the actual data (typically a view)
    panel_size::Tuple{P1, P2} # (rows, columns) in memory
    panel_class::Val{C} # :CM = column major, :CN = row major
    pad_first::Tuple{F1, F2} # # (rows, columns) of padding before the first element of first tile
    pad_last::Tuple{L1, L2} # (rows, columns) of padding after the last element until tile ends

    function Panel{T}(data::D, panel_size::Tuple{P1, P2}, panel_class::Val{C},
            pad_first::Tuple{F1, F2}, pad_last::Tuple{L1, L2}) where {T, D<:Union{Ptr{T}, AbstractArray{T}},
            P1<:Integer, P2<:Integer, C, F1<:Integer, F2<:Integer, L1<:Integer, L2<:Integer}
        new{T,D,P1,P2,C,F1,F2,L1,L2}(data, panel_size, panel_class, pad_first, pad_last)
    end
end

Base.size(x::Panel) = x.panel_size .- x.pad_first .- x.pad_last
Base.eltype(::Type{Panel{T}}) where {T} = T
Base.eltype(::Panel{T}) where {T} = T

"""
    linear_index(panel_class, panel_size, i, j)

Compute a linear index into a panel without padding.
(This function is mainly for internal use. Currently, it does not check bounds,
but this may change.)
"""
@inline linear_index(::Val{:CM}, sz, i, j, pf=static.((0,0))) = i + pf[1] + (j + pf[2] - 1) * sz[1]
@inline linear_index(::Val{:RM}, sz, i, j, pf=static.((0,0))) = linear_index(Val(:CM), reverse(sz), j, i, reverse(pf))

@inline function Base.setindex!(x::Panel, v, i::Integer, j::Integer)
    @boundscheck checkbounds(x, i, j)
    @inbounds x.data[linear_index(x.panel_class, x.panel_size, i, j)] = v
end

@inline function Base.getindex(x::Panel, i::Integer, j::Integer)
    @boundscheck checkbounds(x, i, j)
    @inbounds x.data[linear_index(x.panel_class, x.panel_size, i, j)]
end

const panel_padding_error = ErrorException("Incorrect padding size.")
const panel_index_error = ErrorException("Panel index out of bounds.")

"""
   pie(i, b, p)

Pad-if-equal: returns p if i==b, otherwise 0.
"""
pie(i, b, p) = p === static(0) ? static(0) : i == b ? Int(p) : 0

"""
    panel_view(x, i, j, pad_first, pad_last)

Get a `Panel` view of panel (i,j) in matrix `x`.

If `pad_first` and `pad_last` are not specified, then they are computed based
on the indices (i, j) and the padding of A.
"""
@inline function panel_view(x::PanelMatrix{T}, i::Integer, j::Integer,
        pad_first=(nothing, nothing), pad_last=(nothing, nothing)) where T
    pad_first=something.(pad_first, pie.((i, j), (1, 1), x.pad_first))
    pad_last=something.(pad_last, pie.((i, j), x.n_panels, x.pad_last))
    length = prod(x.panel_size)
    if all(isa.(StaticInteger, typeof.(x.panel_size)))
        length = static(length)
    end
    zeroth = ((i - 1) + (j - 1) * x.panel_stride) * length
    @boundscheck begin
        all(0 .<= pad_first) && all(0 .<= pad_last) && all(pad_first .+ pad_last .<= x.panel_size) || throw(panel_padding_error)
        all(1 .<= (i, j) .<= x.n_panels) || throw(panel_index_error)
    end
    Panel{T}(@inbounds(view(x.data, LengthUnitRange(zeroth, length))), x.panel_size, x.panel_class, pad_first, pad_last)
end

"""
    full_panel_view(x, i, j)

Get a full `Panel` view of panel (i,j) in matrix `x`. A view of an entire panel is
returned even if the matrix does not cover this entire panel.
"""
Base.@propagate_inbounds function full_panel_view(x::PanelMatrix{T}, i::Integer, j::Integer) where T
    panel_view(x, i, j, static.((0, 0)), static.((0, 0)))
end

"""
    get_panel(x, i, j)

Returns a copy of panel(i, j), whith the size depending on the location of the
panel.

    get_panel(x, i, j, pad_first, pad_last)

Returns a copy of panel (i,j) in matrix `x`. The amount of padding is specified
by input arguments, and may be different from the true size of this panel.
If the panel size and all of the padding is statically known, then an `SMatrix`
is returned. Otherwise a normal `Matrix` is returned.
"""
@inline function get_panel(x::PanelMatrix{T}, i::Integer, j::Integer,
        pad_first=(nothing, nothing), pad_last=(nothing, nothing)) where {T}
    v = panel_view(x, i, j, pad_first, pad_last)
    Matrix{T}(v)
end

@inline function get_panel(x::PanelMatrix{T,D,S1,S2,StaticInteger{P1},StaticInteger{P2}},
        i::Integer, j::Integer,
        pad_first::Tuple{StaticInteger{F1}, StaticInteger{F2}},
        pad_last::Tuple{StaticInteger{L1}, StaticInteger{L2}}) where {T,D,S1,S2,P1,P2,F1,F2,L1,L2}
    @boundscheck begin
        all(0 .<= pad_first) && all(0 .<= pad_last) && all(pad_first .+ pad_last .<= x.panel_size) || throw(panel_padding_error)
        all(1 .<= (i, j) .<= x.n_panels) || throw(panel_index_error)
    end
    v = @inbounds panel_view(x, i, j, pad_first, pad_last)
    SMatrix{P1-F1-L1, P2-F2-L2, T, (P1-F1-L1)*(P2-F2-L2)}(v)
end

"""
    get_full_panel(x, i, j)

Get a copy of panel (i, j) in matrix `x`. An entire panel is
returned even if the matrix does not cover this entire panel.
"""
Base.@propagate_inbounds get_full_panel(x, i, j) = get_panel(x, i, j, (static(0), static(0)), (static(0), static(0)))

"""
    set_panel!(x, y, i, j, pad_first, pad_last)

Set panel (i,j) in matrix `x` to `y`. The amount of padding is specified
by input arguments, and may be different from the true size of this panel.
"""
@inline function set_panel!(x::PanelMatrix, y::AbstractMatrix, i::Integer, j::Integer,
        pad_first=(nothing, nothing), pad_last=(nothing, nothing))
    v = panel_view(x, i, j, pad_first, pad_last)
    panel_size = size(v)
    @boundscheck begin
        all(0 .<= pad_first) && all(0 .<= pad_last) && all(pad_first .+ pad_last .<= panel_size) || throw(panel_padding_error)
        checkbounds(y, 1:x.panel_size[1]-pad_first[1]-pad_last[1], x.panel_size[2]-pad_first[2]-pad_last[2])
        all(1 .<= (i, j) .<= x.n_panels) || throw(panel_index_error)
    end
    @inbounds v .= y
end

"""
    set_full_panel!(x, y, i, j)

Set the full panel (i,j) in matrix `x` to `y`. An entire panel is
set even if the matrix does not cover this entire panel.
"""
Base.@propagate_inbounds set_full_panel!(x, y, i, j) = set_panel!(x, y, i, j, (static(0), static(0)), (static(0), static(0)))
