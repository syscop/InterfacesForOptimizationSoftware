"""
A `Panel` represents one of the panels in a Panel matrix.

It is used like a view into the matrix.
"""
struct Panel{T,D,M,N,S,O,P,Q,R} <: AbstractMatrix{T}
    data::D # the actual data (typically a view)
    panel_size::Tuple{M, N} # (rows, columns) in memory
    panel_class::Val{S} # :CM = column major, :CN = row major
    pad_first::Tuple{O, P} # # (rows, columns) of padding before the first element of first tile
    pad_last::Tuple{Q, R} # (rows, columns) of padding after the last element until tile ends

    function Panel{T}(data::D, panel_size::Tuple{M, N}, panel_class::Val{S},
            pad_first::Tuple{O, P}, pad_last::Tuple{Q, R}) where {T, D<:Union{Ptr{T}, AbstractArray{T}},
            M<:Integer, N<:Integer, S, O<:Integer, P<:Integer, Q<:Integer, R<:Integer}
        new{T,D,M,N,S,O,P,Q,R}(data, panel_size, panel_class, pad_first, pad_last)
    end
end

Base.size(x::Panel) = x.panel_size .- x.pad_first .- x.pad_last
Base.eltype(t::Type{Panel{T}}) where {T} = T
Base.eltype(x::Panel{T}) where {T} = T

"""
    linear_index(panel_class, panel_size, i, j)

Compute a linear index into a panel.
(This function is mainly for internal use. Currently, it does not check bounds,
but this may change.)
"""
@inline linear_index(::Val{:CM}, sz, i, j) = i + (j - 1) * sz[1]
@inline linear_index(::Val{:RM}, sz, i, j) = j + (i - 1) * sz[2]

@inline function Base.setindex!(x::Panel, v, i::Integer, j::Integer)
    @boundscheck checkbounds(x, i, j)
    if x.data isa Ptr
        unsafe_store!(x.data, v, linear_index(x.panel_class, x.panel_size, i, j))
    else
        @inbounds x.data[linear_index(x.panel_class, x.panel_size, i, j)] = v
    end
end

@inline function Base.getindex(x::Panel, i::Integer, j::Integer)
    @boundscheck checkbounds(x, i, j)
    if x.data isa Ptr
        unsafe_load(x.data, linear_index(x.panel_class, x.panel_size, i, j))
    else
        @inbounds x.data[linear_index(x.panel_class, x.panel_size, i, j)]
    end
end

"""
    full_panel_view(x, i, j)

Get a `Panel` view of panel (i,j) in matrix `x`. An entire panel is
returned even if the matrix `x` does not cover the entire panel.
"""
@inline function full_panel_view(x::PanelMatrix{T}, i::Integer, j::Integer) where T
    length = prod(x.panel_size)
    if all(isa.(StaticInteger, typeof.(x.panel_size)))
        length = static(length)
    end
    zeroth = ((i - 1) + (j - 1) * x.panel_stride) * length
    Panel{T}(view(x.data, LengthUnitRange(zeroth, length)), x.panel_size, x.panel_class, static.((0, 0)), static.((0, 0)))
end

"""
    unsafe_full_panel_view(x, i, j)

Get an unsafe `Panel` view of panel (i,j) in matrix `x`. An entire panel is
returned even if the matrix `x` does not cover the entire panel.

For performance reasons, a pointer is used instead of a proper object reference.
This means that the original PanelMatrix must be protected from garbage collection
while an unsafe `Panel` is in scope, or memory corruption may occur.

```
   GC.@preserve A let a = unsafe_full_panel_view(A, i, j)
       # do stuff with a, but don't let it escape the scope of the block.
   end
```
"""
@inline function unsafe_full_panel_view(x::PanelMatrix{T}, i::Integer, j::Integer) where T
    length = prod(x.panel_size)
    zeroth = ((i - 1) + (j - 1) * x.panel_stride) * length
    Panel{T}(pointer(x.data, zeroth+1), x.panel_size, x.panel_class, static.((0, 0)), static.((0, 0)))
end
