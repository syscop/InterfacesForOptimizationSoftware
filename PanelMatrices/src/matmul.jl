# Error-throwing helper functions
@noinline errAB(mA, nA, mB, nB, m=true) = throw(DimensionMismatch("$(m ? :matrix : :panel ) A has dimensions ($mA,$nA), $(m ? :matrix : :panel ) B has dimensions ($mB,$nB)"))
@noinline errC(mC, nC, mA, nB, m=true) = throw(DimensionMismatch("$(m ? :matrix : :panel ) C has dimensions $(size(C)), needs ($mA,$nB)"))

function LinearAlgebra.mul!(C::PanelMatrix, A::PanelMatrix, B::PanelMatrix, α::Number=false, β::Number=true)
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

    # For now, we only handle matrices where pad_first is zero.
    # TODO: Handle general case. Check that padding matches.
    @assert A.pad_first == (0, 0)
    @assert B.pad_first == (0, 0)
    @assert C.pad_first == (0, 0)

    # For now, we only handle matrices where pad_last is zero.
    # TODO: Handle general case. Check that padding matches.
    @assert A.pad_last == (0, 0)
    @assert B.pad_last == (0, 0)
    @assert C.pad_last == (0, 0)

    for i in 1 : mA ÷ pmA
        for j in 1 : nB ÷ pnB
            @inbounds begin
                Cp = α*get_full_panel(C, i, j)
                for k = 1 : nA ÷ pnA
                    Cp = Cp .+ β .* (get_full_panel(A, i, k) * get_full_panel(B, k, j))
                end
                set_full_panel!(C, Cp, i, j)
            end
        end
    end
end
