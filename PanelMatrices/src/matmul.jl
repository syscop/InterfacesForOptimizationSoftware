# Error-throwing helper functions
@noinline errAB(mA, nA, mB, nB, m=true) = throw(DimensionMismatch("$(m ? :matrix : :panel ) A has dimensions ($mA,$nA), $(m ? :matrix : :panel ) B has dimensions ($mB,$nB)"))
@noinline errC(mC, nC, mA, nB, m=true) = throw(DimensionMismatch("$(m ? :matrix : :panel ) C has dimensions ($mC,$nC), needs ($mA,$nB)"))

# The SArray multiply in StaticArrays.jl is missing a @fastmath annotation...
@generated function Base.muladd(A::SArray{Tuple{M,K},T,2}, B::SArray{Tuple{K,N},T,2}, C::SArray{Tuple{M,N},T,2}) where {M,N,K,T}
    exprs = [:(+($([:(A[$m,$k]*B[$k,$n]) for k = 1:K]...))) for m = 1:M, n = 1:N]
    return quote
        Base.@_inline_meta
        return $C(@fastmath @inbounds C.data .+ tuple($(exprs...)))
    end
end

# for completeness:
Base.muladd(::SArray{Tuple{M,0},T,2}, ::SArray{Tuple{0,N},T,2}, C::SArray{Tuple{M,N},T,2}) where {M,N,K,T} = C

function LinearAlgebra.mul!(C::PanelMatrix, A::PanelMatrix, B::PanelMatrix, α::Number=static(false))
    mA, nA = size(A)
    mB, nB = size(B)
    mC, nC = size(C)
    if mB != nA
        errAB(mA, nA, mB, nB)
    end
    if mC != mA || nC != nB
        errC(mC, nC, mA, nB)
    end

    # Panel sizes must be compatible.
    # We don't auto-convert because the user should know about the preformance loss.
    pmA, pnA = A.panel_size
    pmB, pnB = B.panel_size
    if pmB != pnA
        errAB(pmA, pnA, pmB, pnB, false)
    end
    pmC, pnC = C.panel_size
    if pmC != pmA || pnC != pnB
        errC(mC, nC, mA, nB, false)
    end

    I, J = C.n_panels
    K = A.n_panels[2]
    @assert K > 1 # TODO: Handle special cases

    pfm, pfn = C.pad_first
    pfk = A.pad_first[2]

    plm, pln = C.pad_last
    plk = A.pad_last[2]

    # TODO: More cache-friendly ordering of panels
    mul_panel!(C, A, B, α, 1, 1, K, pfm, pfn, pfk, static(0), static(0), plk)
    for i = 2 : I - 1
        mul_panel!(C, A, B, α, i, 1, K, static(0), pfn, pfk, static(0), static(0), plk)
    end
    mul_panel!(C, A, B, α, I, 1, K, static(0), pfn, pfk, plm, static(0), plk)
    for j in 2 : J - 1
        mul_panel!(C, A, B, α, 1, j, K, pfm, static(0), pfk, static(0), static(0), plk)
        for i = 2 : I - 1
            mul_panel!(C, A, B, α, i, j, K, static(0), static(0), pfk, static(0), static(0), plk)
        end
        mul_panel!(C, A, B, α, I, j, K, static(0), static(0), pfk, plm, static(0), plk)
    end
    mul_panel!(C, A, B, α, 1, J, K, pfm, static(0), pfk, static(0), pln, plk)
    for i = 2 : I - 1
        mul_panel!(C, A, B, α, i, J, K, static(0), static(0), pfk, static(0), pln, plk)
    end
    mul_panel!(C, A, B, α, I, J, K, static(0), static(0), pfk, plm, pln, plk)

    return C
end

"""
Compute panel (i,j) of C = A*B.

The amount of padding is passed explicitly as arguments, rather than being
computed as a function of (i,j)
"""
function mul_panel!(C, A, B, α, i, j, K, pfm, pfn, pfk, plm, pln, plk)
    Cp = α .* @inbounds(get_panel(C, i, j, (pfm, pfn), (plm, pln)))
    Aik = @inbounds get_panel(A, i, 1, (pfm, pfk), (plm, static(0)))
    Bkj = @inbounds get_panel(B, 1, j, (pfk, pfn), (static(0), pln))
    Cp = muladd(Aik, Bkj, Cp)
    for k = 2 : K - 1
        Aik = @inbounds get_panel(A, i, k, (pfm, static(0)), (plm, static(0)))
        Bkj = @inbounds get_panel(B, k, j, (static(0), pfn), (static(0), pln))
        Cp = muladd(Aik, Bkj, Cp)
    end
    Aik = @inbounds get_panel(A, i, K, (pfm, static(0)), (plm, plk))
    Bkj = @inbounds get_panel(B, K, j, (static(0), pfn), (plk, pln))
    Cp = muladd(Aik, Bkj, Cp)
    @inbounds set_panel!(C, Cp, i, j, (pfm, pfn), (plm, pln))
end
