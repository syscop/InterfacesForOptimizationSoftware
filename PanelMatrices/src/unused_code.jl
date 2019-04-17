"""
    unsafe_full_panel_view(x, i, j)

NOTE: This function might be removed in the future. Use `full_panel_view` where possible.

Get an unsafe `Panel` view of panel (i,j) in matrix `x`. An entire panel is
returned even if the matrix does not cover this entire panel.

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
    @boundscheck all(1 .<= (i, j) .* x.panel_size .<= x.size .+ x.pad_first .+ x.pad_last) || throw(panel_index_error)
    Panel{T}(pointer(x.data, zeroth+1), x.panel_size, x.panel_class, static.((0, 0)), static.((0, 0)))
end
