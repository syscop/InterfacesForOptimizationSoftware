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
        # panel from unsafe_full_panel_view
        unsafe_store!(x.data, v, linear_index(x.panel_class, x.panel_size, i, j))
    else
        @inbounds x.data[linear_index(x.panel_class, x.panel_size, i, j)] = v
    end
end

@inline function Base.getindex(x::Panel, i::Integer, j::Integer)
    @boundscheck checkbounds(x, i, j)
    if x.data isa Ptr
        # panel from unsafe_full_panel_view
        unsafe_load(x.data, linear_index(x.panel_class, x.panel_size, i, j))
    else
        @inbounds x.data[linear_index(x.panel_class, x.panel_size, i, j)]
    end
end

const panel_index_error = ErrorException("Panel index out of bounds.")

"""
    full_panel_view(x, i, j)

Get a `Panel` view of panel (i,j) in matrix `x`. An entire panel is
returned even if the matrix does not cover this entire panel.
"""
@inline function full_panel_view(x::PanelMatrix{T}, i::Integer, j::Integer) where T
    length = prod(x.panel_size)
    if all(isa.(StaticInteger, typeof.(x.panel_size)))
        length = static(length)
    end
    zeroth = ((i - 1) + (j - 1) * x.panel_stride) * length
    @boundscheck all(1 .<= (i, j) .* x.panel_size .<= x.size .+ x.pad_first .+ x.pad_last) || throw(panel_index_error)
    Panel{T}(@inbounds(view(x.data, LengthUnitRange(zeroth, length))), x.panel_size, x.panel_class, static.((0, 0)), static.((0, 0)))
end

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

"""
    get_panel(x, i, j, pad_top, pad_bottom, pad_left, pad_right)

Returns a copy of panel (i,j) in matrix `x`. The amount of padding is specified
by input arguments, and may be different from the true size of this panel.
If the panel size and all of the padding is statically known, then an `SMatrix`
is returned. Otherwise a normal `Matrix` is returned.
"""
@generated function get_panel(x::PanelMatrix{T,D,U,V,StaticInteger{M},StaticInteger{N}},
        i::Integer, j::Integer,
        ::StaticInteger{PT}, ::StaticInteger{PB},
        ::StaticInteger{PL}, ::StaticInteger{PR}) where {T,D,U,V,M,N,PT,PB,PL,PR}
    quote
        Base.@_inline_meta
        @boundscheck all(0 .<= (PT, PB, PL, PR)) && all((PT, PL) .+ (PB, PR) .<= (M,N)) || throw(panel_index_error)
        v = full_panel_view(x, i, j)
        SMatrix{$(M-PT-PB),$(N-PL-PR),$T,$((M-PT-PB)*(N-PL-PR))}(
            @inbounds $(Expr(:tuple, vec([:(v[$i,$j]) for i=PT+1:M-PB, j=PL+1:N-PR])...)))
    end
end
function get_panel(x::PanelMatrix, i::Integer, j::Integer, pt::Integer, pb::Integer, pl::Integer, pr::Integer)
    v = full_panel_view(x, i, j)
    (m,n) = size(v)
    @boundscheck all(0 .<= (pt, pb, pl, pr)) && all((pt, pl) .+ (pb, pr) .<= (m, n)) || throw(panel_index_error)
    @inbounds v[pt+1:m-pb, pl+1:n-pr]
end

"""
    get_panel(x, i, j)

Returns a copy of panel(i, j), whith the size depending on the location of the
panel.
"""
# TODO

"""
    get_full_panel(x, i, j)

Get a copy of panel (i, j) in matrix `x`. An entire panel is
returned even if the matrix does not cover this entire panel.
"""
get_full_panel(x, i, j) = get_panel(x, i, j, static(0), static(0), static(0), static(0))

"""
    set_full_panel!(x, y, i, j)

Get a copy of panel (i,j) in matrix `x`. An entire panel is
returned even if the matrix does not cover this entire panel.
"""
@generated function set_full_panel!(x::PanelMatrix{T,D,U,V,StaticInteger{M},StaticInteger{N}},
        y::AbstractArray, i::Integer, j::Integer) where {T,D,U,V,M,N}
    quote
        Base.@_inline_meta
        @boundscheck checkbounds(y, 1:M, 1:N)
        v = full_panel_view(x, i, j)
        @inbounds ($(Expr(:tuple, vec([:(v[$i,$j]) for i=1:M, j=1:N])...))) = ($(Expr(:tuple, vec([:(y[$i,$j]) for i=1:M, j=1:N])...)))
    end
end
