module PanelMatrices

using StaticNumbers
using StaticArrays
#using UnsafeArrays # maybe we don't need this, but just use the same idea.

export PanelMatrix

const panel_stride_error = ErrorException("Panel stride is less than minimum required.")

# TODO: Check if Julia specializes on the tuples, or if they need to be taken apart?

"""
A `PanelMatrix` is a panel-major matrix, or a view into a panel-major matrix.

It contains a `.data` field which holds the actual data. This must be an array
that maps to continuous memory storage. (For efficiency, panels may be accessed
using hand-coded assembly that makes this assumtion.)

Care should be taken if mutating the `.data` array. Changing its size
(e.g. using `push!` or `pop!`) may lead to data corruption.

Memory is allocated to cover the entire panel, even in cases where only parts
of a panel is covered by the matrix.
"""
struct PanelMatrix{T,D,U,V,M,N,S,O,P,Q,R} <: AbstractMatrix{T}
    data::D # The actual data <: AbstractVector{T}
    size::Tuple{U, V} # (rows, columns) in the matrix
    panel_size::Tuple{M, N} # (rows, columns) in each panel
    panel_class::Val{S}  # :CM = column major, :RM = row major
    pad_first::Tuple{O, P} # (rows, columns) of padding before the first element of first tile
    pad_last::Tuple{Q, R} # (rows, columns) of padding after the last element until last tile ends
    panel_stride::Int # number of panels per column in data

    # The inner constructor checks the size of the data vector.
    # Use @inbounds to bypass size and stride check.
    # (Mutating the size of .data may lead to memory corroption. Don't do that!)
    # pad_last is computed from pad_first, size and panel_size. Only its type can be specified.
    function PanelMatrix{T,D}(data::D, size::Tuple{U, V}, panel_size::Tuple{M, N},
           panel_class::Val{S}, pad_first::Tuple{O, P}, pad_last_types::Tuple{<:Type, <:Type},
           panel_stride::Int) where {T, D<:AbstractVector{T}, U<:Integer, V<:Integer, M<:Integer, N<:Integer, S, O<:Integer, P<:Integer}
        @boundscheck begin
            panel_stride * panel_size[1] >= size[1] + pad_first[1] || throw(panel_stride_error)
            n_panels = cld.(pad_first .+ size, panel_size)
            last_panel = n_panels[1] + (n_panels[2]-1)*panel_stride
            checkbounds(data, Base.OneTo(last_panel*prod(panel_size)))
        end
        pad_last = convert.(pad_last_types, mod.(panel_size .- pad_first .- size, panel_size))
        (Q, R) = typeof.(pad_last)
        new{T,D,U,V,M,N,S,O,P,Q,R}(data, size, panel_size, panel_class, pad_first, pad_last, panel_stride)
    end
end

# TODO: We could use keyword args, but first benchmark without and check that
#       there's no performance issue.

"""
    PanelMatrix{T}(undef, size, panel_size, panel_class, pad_first, pad_last_types)

Create an uninitialized `PanelMatrix` for element type `T`, with the given `size`.
"""
function PanelMatrix{T}(
       ::UndefInitializer,
       size::Tuple{<:Integer, <:Integer},
       panel_size::Tuple{<:Integer, <:Integer} = default_panelsize(eltype(data)),
       panel_class::Val = default_panelclass(T),
       pad_first::Tuple{<:Integer, <:Integer} = static.((0, 0)),
       pad_last_types::Tuple{<:Type, <:Type} = map(t isa StaticInteger ? StaticInteger : Int, typeof(size))
       ) where {T}
    panel_stride = cld(size[1]+pad_first[1], panel_size[1])
    n_panels = cld.(pad_first .+ size, panel_size)
    last_panel = n_panels[1] + (n_panels[2]-1)*panel_stride
    data = Vector{T}(undef, last_panel*prod(panel_size))
    PanelMatrix{T, Vector{T}}(data, size, panel_size, panel_class, pad_first, pad_last_types, panel_stride)
end

"""
    PanelMatrix(x, panel_size, panel_class, pad_first, pad_last_types)

Create a `PanelMatrix` from the matrix `x` (making a copy).
"""
function PanelMatrix(x::AbstractMatrix,
        panel_size::Tuple{<:Integer, <:Integer} = default_panelsize(eltype(x)),
        panel_class::Val = default_panelclass(T),
        pad_first::Tuple{<:Integer, <:Integer} = static.((0, 0)),
        pad_last_types::Tuple{<:Type, <:Type} = (Int, Int)
        )
    z = PanelMatrix{eltype(x)}(undef, size(x), panel_size, panel_class, pad_first, pad_last_types)
    #  TODO: implement z .= x
    #  TODO: fill padding with zeros ?
    for i in eachindex(x)
        z[i] = x[i]
    end
    return z
end

default_panelsize(T) = static.((4, 4))

default_panelclass(T) = Val(:CM)

include("Panels.jl")

@inline function Base.getindex(x::PanelMatrix, i, j)
    @boundscheck checkbounds(x, i, j)
    (p, q) = divrem.((i, j) .+ x.pad_first .- 1, x.panel_size)
    @inbounds full_panel_view(x, p[1]+1, q[1]+1)[p[2]+1, q[2]+1]
end

@inline function Base.setindex!(x::PanelMatrix, v, i, j)
    @boundscheck checkbounds(x, i, j)
    (p, q) = divrem.((i, j) .+ x.pad_first .- 1, x.panel_size)
    @inbounds full_panel_view(x, p[1]+1, q[1]+1)[p[2]+1, q[2]+1] = v
end

Base.size(x::PanelMatrix) = x.size
Base.eltype(x::PanelMatrix) = eltype(x.data)

# TODO: Base.similar

end # module
