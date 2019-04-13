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

struct UnsafePanel{T,M,N,S,O,P,Q,R} <: AbstractMatrix{T}
    data::Ptr{T} # pointer to the actual data (L elements of T)

    # These will always be size zero.
    panelsize::Tuple{M, N} # (rows, columns) in memory
    panel_class::Val{S} # :CM = column major, :CN = row major

    # These may be Int or StaticInteger, depending on tradeoff between compile and runtime.
    pad_first::Tuple{O, P} # # (rows, columns) of padding before the first element of first tile
    pad_last::Tuple{Q, R} # (rows, columns) of padding after the last element until tile ends

    function UnsafePanel(data::Ptr{T}, panelsize::Tuple{M, N}, panel_class::Val{S},
            pad_first::Tuple{O, P}, pad_last::Tuple{Q, R}) where {T,
            M<:Integer, N<:Integer, S, O<:Integer, P<:Integer, Q<:Integer, R<:Integer}
        new{T,M,N,S,O,P,Q,R}(data, panelsize, panel_class, pad_first, pad_last)
    end
end

Base.size(x::UnsafePanel) = x.panelsize .- x.pad_first .- x.pad_last
Base.eltype(t::Type{UnsafePanel{T}}) where {T} = T
Base.eltype(x::UnsafePanel{T}) where {T} = T

# TODO: Check panel_class
function Base.setindex!(x::UnsafePanel{T}, v::T, i::Integer, j::Integer) where T
     @boundscheck checkbounds(x, i, j)
     unsafe_store!(x.data, v, i + (j - 1) * x.panelsize[1])
end

# TODO: Check panel_class
function Base.getindex(x::UnsafePanel, i::Integer, j::Integer)
    @boundscheck checkbounds(x, i, j)
    unsafe_load(x.data, i + (j - 1) * x.panelsize[1])
end

function unsafe_get_full_panel(x::PanelMatrix, i::Integer, j::Integer)
    offset = ((i - 1) + (j - 1) * x.panel_stride) * prod(x.panel_size) + 1
    UnsafePanel(pointer(x.data, offset), x.panel_size, x.panel_class, static.((0, 0)), static.((0, 0)))
end
