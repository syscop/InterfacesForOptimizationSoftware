module PanelMatrices

using StaticNumbers
using StaticArrays
#using UnsafeArrays # maybe we don't need this, but just use the same idea.

export PanelMatrix

const panel_stride_error = ErrorException("Panel stride is less than minimum required.")
const too_small_error = ErrorException("PanelMatrix size must either be static or greater than the panel size.")

# TODO: Move maybe_static() to StaticNumbers and make a macro for it.
"""
    maybe_static(f, args...)

Returns `static(f(args...))`, if all of args are `Static`. If any of the args
is not `Static`, then `f(args...)` is returned unchanged.
"""
@inline maybe_static(f::F, args...) where {F} = f(args...)
@inline maybe_static(f::F, args::Static...) where {F} = static(f(args...))

# The inner PanelMatrix constructor checks the size of the data vector.
# Use @inbounds to bypass size and stride check.
# (Mutating the size of .data may lead to memory corroption. Don't do that!)
# We force the size to be static when smaller than or equal to one panel.
# All subsequent code can assume that if size[i] isa Int then that Int
# is not smaller than panel_size[i] and therefore that pat_first[i] and
# pad_last[i] never apply to the same panel.
# pad_last is computed from pad_first, size and panel_size. Only its type can be specified.
# n_panels is computed from pad_first, size and panel_size. It will be Static if all of those are.

"""
A `PanelMatrix` is a panel-major matrix, or a view into a panel-major matrix.

The size of a PanelMatrix can either be Int or StaticInteger. In either case,
the  matrix cannot be resized. (Using Int requires fewer methods to be compiled
when matrices of varying sizes are used.)

If either the width or height of a PanelMatrix is smaller than the default panel
size, then it is recommended to use a StaticInteger for that element of the size,
and to reduce the panel size so that it does not exceed the matrix size.

A PanelMatrix contains a `.data` field which holds the actual data. This must
be an array that maps to continuous memory storage. (For efficiency, panels may
be accessed using hand-coded assembly that makes this assumtion.)

Care should be taken if mutating the `.data` array. Changing its size
(e.g. using `push!` or `pop!`) may lead to data corruption.

Memory is allocated to cover the entire panel, even in cases where only parts
of a panel is covered by the matrix.
"""
struct PanelMatrix{T,D,S1,S2,P1,P2,C,F1,F2,L1,L2,R,N1,N2} <: AbstractMatrix{T}
    data::D # The actual data <: AbstractVector{T}
    size::Tuple{S1, S2} # (rows, columns) in the matrix
    panel_size::Tuple{P1, P2} # (rows, columns) in each panel
    panel_class::Val{C}  # :CM = column major, :RM = row major
    pad_first::Tuple{F1, F2} # (rows, columns) of padding before the first element of first tile
    pad_last::Tuple{L1, L2} # (rows, columns) of padding after the last element until last tile ends
    panel_stride::R # number of panels per column in data
    n_panels::Tuple{N1, N2} # (rows, columns) of panels in the matrix

    function PanelMatrix{T,D}(data::D, size::Tuple{S1, S2}, panel_size::Tuple{P1, P2},
           panel_class::Val{C}, pad_first::Tuple{F1, F2}, pad_last_types::Tuple{<:Type, <:Type},
           panel_stride::R) where {T, D<:AbstractVector{T}, S1<:Integer, S2<:Integer, P1<:Integer, P2<:Integer, C, F1<:Integer, F2<:Integer, R<:Integer}

        n_panels = maybe_static.((f,s,p) -> cld.(f+s,p), pad_first, size, panel_size)
        @boundscheck begin
            panel_stride * panel_size[1] >= size[1] + pad_first[1] || throw(panel_stride_error)
            last_panel = n_panels[1] + (n_panels[2]-1)*panel_stride
            checkbounds(data, Base.OneTo(last_panel*prod(panel_size)))
        end
        any(.!(isa.(StaticInteger, typeof.(size))) .& (size .<= panel_size)) && throw(too_small_error)
        pad_last = convert.(pad_last_types, panel_size .* n_panels .- pad_first .- size)
        (L1, L2) = typeof.(pad_last)
        (N1, N2) = typeof.(n_panels)
        new{T,D,S1,S2,P1,P2,C,F1,F2,L1,L2,R,N1,N2}(data, size, panel_size, panel_class, pad_first, pad_last, panel_stride, n_panels)
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
       panel_size::Tuple{<:Integer, <:Integer} = static.(min.(size, default_panelsize(T))),
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
        panel_size::Tuple{<:Integer, <:Integer} = static.(min.(size(x), default_panelsize(eltype(x)))),
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
