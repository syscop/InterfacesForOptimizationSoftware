# PanelMatrices.UnsafePanel is unexported and should be used with care.
"""
An `UnsafePanel` is similar to a view into a column major matrix.

For performance reasons, a pointer is used instead of a proper object reference.
This means that the original PanelMatrix must be protected from garbage collection
while an UnsafePanel is in scope, or memory corruption may occur.

```
   GC.@preserve A let a = unsafe_get_panel(A, i, j)
       # do stuff with a, but don't let it escape the scope of the block.
   end
```
"""

struct UnsafePanel{T,U,V,O,P,S,M,N,L} <: AbstractMatrix{T}
    data::Ptr{T} # pointer to the actual data (L elements of T)

    # These may be Int or StaticInteger, depending on tradeoff between compile and runtime.
    sub_first::Tuple{O, P} # (row, column) of first element, within the tile
    sub_last::Tuple{Q, R} # (row, column) of last element, within the tile

    # These will always be size zero.
    panel_class::Val{S} # :CM = column major, :CN = row major
    fullsize::Tuple{StaticInteger{M}, StaticInteger{N}} # (rows, columns) in memory
    datalength::StaticInteger{L} # prod(fullsize)
end

Base.size(t::Type{UnsafePanel}) = t.size
Base.size(x::UnsafePanel) = x.size
Base.eltype(t::Type{UnsafePanel{T}}) where {T} = T
Base.eltype(x::UnsafePanel{T}) where {T} = T

function Base.setindex!(x::UnsafePanel{T}, v::T, i::Integer, j::Integer) where T
     @boundscheck checkbounds(x, i, j)
     unsafe_store!(x.data, v, i + (j - 1) * x.fullsize[1])
end

function Base.getindex(x::UnsafePanel, i::Integer, j::Integer)
    @boundscheck checkbounds(x, i, j)
    unsafe_load(x.data, i + (j - 1) * x.fullsize[1])
end
