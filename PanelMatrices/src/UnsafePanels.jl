# PanelMatrices.UnsafePanel is unexported and should be used with care.
"""
An `UnsafePanel` is similar to a view into a panel matrix. It is intended for
internal use by PanelMatrices.

For performance reasons, a pointer is used instead of a proper object reference.
This means that the original PanelMatrix must be protected from garbage collection
while an UnsafePanel is in scope, or memory corruption may occur.

```
   GC.@preserve A let a = unsafe_get_panel(A, i, j)
       # do stuff with a, but don't let it escape the scope of the block.
   end
```
"""

struct UnsafePanel{T,M,N,S,O,P,Q,R} <: AbstractMatrix{T}
    data::Ptr{T} # pointer to the actual data (L elements of T)

    # These will always be size zero.
    panel_size::Tuple{M, N} # (rows, columns) in memory
    panel_class::Val{S} # :CM = column major, :CN = row major

    # These may be Int or StaticInteger, depending on tradeoff between compile and runtime.
    pad_first::Tuple{O, P} # # (rows, columns) of padding before the first element of first tile
    pad_last::Tuple{Q, R} # (rows, columns) of padding after the last element until tile ends

    function UnsafePanel(data::Ptr{T}, panel_size::Tuple{M, N}, panel_class::Val{S},
            pad_first::Tuple{O, P}, pad_last::Tuple{Q, R}) where {T,
            M<:Integer, N<:Integer, S, O<:Integer, P<:Integer, Q<:Integer, R<:Integer}
        new{T,M,N,S,O,P,Q,R}(data, panel_size, panel_class, pad_first, pad_last)
    end
end

Base.size(x::UnsafePanel) = x.panel_size .- x.pad_first .- x.pad_last
Base.eltype(t::Type{UnsafePanel{T}}) where {T} = T
Base.eltype(x::UnsafePanel{T}) where {T} = T

"""
    linear_index(panel_class, panel_size, i, j)

Compute a linear index into a panel.
(This function is mainly for internal use. Currently, it does not check bounds,
but this may change.)
"""
@inline linear_index(::Val{:CM}, sz, i, j) = i + (j - 1) * sz[1]
@inline linear_index(::Val{:RM}, sz, i, j) = j + (i - 1) * sz[2]

function Base.setindex!(x::UnsafePanel{T}, v::T, i::Integer, j::Integer) where T
     @boundscheck checkbounds(x, i, j)
     unsafe_store!(x.data, v, linear_index(x.panel_class, x.panel_size, i, j))
end

function Base.getindex(x::UnsafePanel, i::Integer, j::Integer)
    @boundscheck checkbounds(x, i, j)
    unsafe_load(x.data, linear_index(x.panel_class, x.panel_size, i, j))
end

"""
    unsafe_full_panel_view(x, i, j)

Get an `UnsafePanel` view of panel (i,j) in matrix `x`. An entire panel is
returned even if the matrix `x` does not cover the entire panel.
"""
function unsafe_full_panel_view(x::PanelMatrix, i::Integer, j::Integer)
    offset = ((i - 1) + (j - 1) * x.panel_stride) * prod(x.panel_size) + 1
    UnsafePanel(pointer(x.data, offset), x.panel_size, x.panel_class, static.((0, 0)), static.((0, 0)))
end
