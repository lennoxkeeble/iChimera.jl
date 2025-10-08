module MultipoleDerivs
using ChimeraMultipoles
using StaticArrays
using ..SymmetricTensors
using ..CoordinateDerivs
using ..HarmonicCoordDerivs

# define mass-ratio parameter
@inline η(q::Float64)::Float64 = q / ((1+q)^2)   # q = SCO / MBH (mass ratio < 1)
@inline Mass_quad_prefactor(q::Float64)::Float64 = η(q) * (1.0 + q)
@inline Mass_oct_prefactor(q::Float64)::Float64 = η(q) * (1.0 - q) # there should be an overall minus sign ahead of Mijk, this minus sign is already absorbed in the expressions below
@inline Mass_hex_prefactor(q::Float64)::Float64 = η(q) * (1.0 + q)
@inline Current_quad_prefactor(q::Float64)::Float64 = η(q) * (1.0 - q) # there should be an overall minus sign ahead of Mijk, this minus sign is already absorbed in the expressions below
@inline Current_oct_prefactor(q::Float64)::Float64 = η(q) * (1.0 + q)

function compute_SF_moments!(q::Float64, Mij5::MMatrix{3, 3, Float64, 9}, Mij6::MMatrix{3, 3, Float64, 9}, Mij7::MMatrix{3, 3, Float64, 9}, Mij8::MMatrix{3, 3, Float64, 9}, dxmMij5::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, dxmMij6::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, dxmMij7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijk7::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijk8::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, dxmMijk7::MArray{Tuple{3, 3, 3, 3}, Float64, 4, 81}, Sij5::MMatrix{3, 3, Float64, 9}, Sij6::MMatrix{3, 3, Float64, 9}, dxmSij5::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64}, d5x::MVector{3, Float64}, d6x::MVector{3, Float64}, d7x::MVector{3, Float64}, d8x::MVector{3, Float64}, d9x::MVector{3, Float64}, OnePN::Float64, TwoPN::Float64, TwoPointFivePN::Float64)
    if OnePN != 0.0 && OnePN != 1.0
        error("OnePN must be either 0.0 or 1.0")
    end

    if TwoPN != 0.0 && TwoPN != 1.0
        error("TwoPN must be either 0.0 or 1.0")
    end

    if TwoPointFivePN != 0.0 && TwoPointFivePN != 1.0
        error("TwoPointFivePN must be either 0.0 or 1.0")
    end

    # compute necessary auxiliary quantities
    m = 1.0 + q;
    ν = q / ((1+q)^2);

    # ChimeraMultipoles.SijDerivs.Sij2(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
    # ChimeraMultipoles.MijkDerivs.Mijk3(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
    # ChimeraMultipoles.MijDerivs.Mij2(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
    

    @inbounds Threads.@threads for (i,j) in collect(Iterators.product(1:3, 1:3))
        Mij5[i,j] = Mass_quad_prefactor(q) * ChimeraMultipoles.MijDerivs.Mij5(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        Mij6[i,j] = Mass_quad_prefactor(q) * ChimeraMultipoles.MijDerivs.Mij6(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        Mij7[i,j] = Mass_quad_prefactor(q) * ChimeraMultipoles.MijDerivs.Mij7(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        Mij8[i,j] = Mass_quad_prefactor(q) * ChimeraMultipoles.MijDerivs.Mij8(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
        
        Sij5[i,j] = Current_quad_prefactor(q) * ChimeraMultipoles.SijDerivs.Sij5(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
        Sij6[i,j] = Current_quad_prefactor(q) * ChimeraMultipoles.SijDerivs.Sij6(i, j, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)


        for k in 1:3
            Mijk7[i,j,k] = Mass_oct_prefactor(q) * ChimeraMultipoles.MijkDerivs.Mijk7(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
            Mijk8[i,j,k] = Mass_oct_prefactor(q) * ChimeraMultipoles.MijkDerivs.Mijk8(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)

            dxmMij5[i,j,k] = Mass_quad_prefactor(q) * ChimeraMultipoles.dxMijDerivs.dxkMij5(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
            dxmMij6[i,j,k] = Mass_quad_prefactor(q) * ChimeraMultipoles.dxMijDerivs.dxkMij6(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)
            dxmMij7[i,j,k] = Mass_quad_prefactor(q) * ChimeraMultipoles.dxMijDerivs.dxkMij7(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN, m, ν)


            dxmSij5[i,j,k] = Current_quad_prefactor(q) * ChimeraMultipoles.dxSijDerivs.dxkSij5(i, j, k, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)

            for l = 1:3
                dxmMijk7[i,j,k,l] = Mass_oct_prefactor(q) * ChimeraMultipoles.dxMijkDerivs.dxlMijk7(i, j, k, l, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, m, ν)
            end
        end
    end
end

function compute_WF_moments!(q::Float64, Mij2::MMatrix{3, 3, Float64, 9}, Mijk3::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, Mijkl4::MArray{Tuple{3, 3, 3, 3}, Float64, 4, 81}, Sij2::MMatrix{3, 3, Float64, 9}, Sijk3::MArray{Tuple{3, 3, 3}, Float64, 3, 27}, x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})
    # fill independent components
    Mij2[1, 2] = MultipoleDerivs.D2_M12(x, dx, d2x)
    Mij2[1, 3] = MultipoleDerivs.D2_M13(x, dx, d2x)
    Mij2[2, 3] = MultipoleDerivs.D2_M23(x, dx, d2x)
    Mij2[1, 1] = MultipoleDerivs.D2_M11(x, dx, d2x)
    Mij2[2, 2] = MultipoleDerivs.D2_M22(x, dx, d2x)
    Mij2[3, 3] = MultipoleDerivs.D2_M33(x, dx, d2x)

    Mijk3[1, 1, 1] = MultipoleDerivs.D3_M111(x, dx, d2x, d3x)
    Mijk3[1, 1, 2] = MultipoleDerivs.D3_M112(x, dx, d2x, d3x)
    Mijk3[1, 2, 2] = MultipoleDerivs.D3_M122(x, dx, d2x, d3x)
    Mijk3[1, 1, 3] = MultipoleDerivs.D3_M113(x, dx, d2x, d3x)
    Mijk3[1, 3, 3] = MultipoleDerivs.D3_M133(x, dx, d2x, d3x)
    Mijk3[1, 2, 3] = MultipoleDerivs.D3_M123(x, dx, d2x, d3x)
    Mijk3[2, 2, 2] = MultipoleDerivs.D3_M222(x, dx, d2x, d3x)
    Mijk3[2, 2, 3] = MultipoleDerivs.D3_M223(x, dx, d2x, d3x)
    Mijk3[2, 3, 3] = MultipoleDerivs.D3_M233(x, dx, d2x, d3x)
    Mijk3[3, 3, 3] = MultipoleDerivs.D3_M333(x, dx, d2x, d3x)

    Mijkl4[1, 1, 1, 1] = MultipoleDerivs.D4_M1111(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 1, 1, 2] = MultipoleDerivs.D4_M1112(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 1, 2, 2] = MultipoleDerivs.D4_M1122(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 2, 2, 2] = MultipoleDerivs.D4_M1222(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 1, 1, 3] = MultipoleDerivs.D4_M1113(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 1, 3, 3] = MultipoleDerivs.D4_M1133(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 3, 3, 3] = MultipoleDerivs.D4_M1333(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 1, 2, 3] = MultipoleDerivs.D4_M1123(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 2, 2, 3] = MultipoleDerivs.D4_M1223(x, dx, d2x, d3x, d4x)
    Mijkl4[1, 2, 3, 3] = MultipoleDerivs.D4_M1233(x, dx, d2x, d3x, d4x)
    Mijkl4[2, 2, 2, 2] = MultipoleDerivs.D4_M2222(x, dx, d2x, d3x, d4x)
    Mijkl4[2, 2, 2, 3] = MultipoleDerivs.D4_M2223(x, dx, d2x, d3x, d4x)
    Mijkl4[2, 2, 3, 3] = MultipoleDerivs.D4_M2233(x, dx, d2x, d3x, d4x)
    Mijkl4[2, 3, 3, 3] = MultipoleDerivs.D4_M2333(x, dx, d2x, d3x, d4x)
    Mijkl4[3, 3, 3, 3] = MultipoleDerivs.D4_M3333(x, dx, d2x, d3x, d4x)

    Sij2[1, 2] = MultipoleDerivs.D2_S12(x, dx, d2x, d3x)
    Sij2[1, 3] = MultipoleDerivs.D2_S13(x, dx, d2x, d3x)
    Sij2[2, 3] = MultipoleDerivs.D2_S23(x, dx, d2x, d3x)
    Sij2[1, 1] = MultipoleDerivs.D2_S11(x, dx, d2x, d3x)
    Sij2[2, 2] = MultipoleDerivs.D2_S22(x, dx, d2x, d3x)
    Sij2[3, 3] = MultipoleDerivs.D2_S33(x, dx, d2x, d3x)

    Sijk3[1, 1, 1] = MultipoleDerivs.D3_S111(x, dx, d2x, d3x, d4x)
    Sijk3[1, 1, 2] = MultipoleDerivs.D3_S112(x, dx, d2x, d3x, d4x)
    Sijk3[1, 2, 2] = MultipoleDerivs.D3_S122(x, dx, d2x, d3x, d4x)
    Sijk3[1, 1, 3] = MultipoleDerivs.D3_S113(x, dx, d2x, d3x, d4x)
    Sijk3[1, 3, 3] = MultipoleDerivs.D3_S133(x, dx, d2x, d3x, d4x)
    Sijk3[1, 2, 3] = MultipoleDerivs.D3_S123(x, dx, d2x, d3x, d4x)
    Sijk3[2, 2, 2] = MultipoleDerivs.D3_S222(x, dx, d2x, d3x, d4x)
    Sijk3[2, 2, 3] = MultipoleDerivs.D3_S223(x, dx, d2x, d3x, d4x)
    Sijk3[2, 3, 3] = MultipoleDerivs.D3_S233(x, dx, d2x, d3x, d4x)
    Sijk3[3, 3, 3] = MultipoleDerivs.D3_S333(x, dx, d2x, d3x, d4x)

    # having filled the independent components, now fill entire matrices
    SymmetricTensors.SymmetrizeTwoIndexTensor!(Mij2)
    SymmetricTensors.SymmetrizeThreeIndexTensor!(Mijk3)
    SymmetricTensors.SymmetrizeFourIndexTensor!(Mijkl4)
    SymmetricTensors.SymmetrizeTwoIndexTensor!(Sij2)
    SymmetricTensors.SymmetrizeThreeIndexTensor!(Sijk3)

    Mij2 .= MultipoleDerivs.Mass_quad_prefactor(q) * Mij2
    Mijk3 .= MultipoleDerivs.Mass_oct_prefactor(q) * Mijk3
    Mijkl4 .= MultipoleDerivs.Mass_hex_prefactor(q) * Mijkl4
    Sij2 .= MultipoleDerivs.Current_quad_prefactor(q) * Sij2
    Sijk3 .= MultipoleDerivs.Current_oct_prefactor(q) * Sijk3;
end

@inline D2_M12(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64})::Float64 = 2*dx[1]*dx[2] + x[2]*d2x[1] + x[1]*d2x[2]

@inline D2_M13(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64})::Float64 = 2*dx[1]*dx[3] + x[3]*d2x[1] + x[1]*d2x[3]

@inline D2_M23(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64})::Float64 = 2*dx[2]*dx[3] + x[3]*d2x[2] + x[2]*d2x[3]

@inline D2_M11(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64})::Float64 = (-2*(-2*dx[1]^2 + dx[2]^2 + dx[3]^2 - 2*x[1]*d2x[1] + x[2]*d2x[2] + x[3]*d2x[3]))/3.

@inline D2_M22(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64})::Float64 = (-2*(dx[1]^2 - 2*dx[2]^2 + dx[3]^2 + x[1]*d2x[1] - 2*x[2]*d2x[2] + x[3]*d2x[3]))/3.

@inline D2_M33(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64})::Float64 = (-2*(dx[1]^2 + dx[2]^2 - 2*dx[3]^2 + x[1]*d2x[1] + x[2]*d2x[2] - 2*x[3]*d2x[3]))/3.

@inline D3_M111(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (3*(-4*dx[1]^3 + 6*x[3]*dx[3]*d2x[1] + 6*x[1]*dx[2]*d2x[2] + 6*x[1]*dx[3]*d2x[3] + 6*dx[1]*(dx[2]^2 + dx[3]^2 - 2*x[1]*d2x[1] + x[2]*d2x[2] + x[3]*d2x[3]) - 2*x[1]^2*d3x[1] + x[2]^2*d3x[1] + x[3]^2*d3x[1] + 2*x[2]*(3*dx[2]*d2x[1] + x[1]*d3x[2]) + 2*x[1]*x[3]*d3x[3]))/5.

@inline D3_M112(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (-24*dx[1]^2*dx[2] + 6*dx[2]^3 + 6*x[3]*dx[3]*d2x[2] - 24*dx[1]*(x[2]*d2x[1] + x[1]*d2x[2]) + 6*x[2]*dx[3]*d2x[3] + 6*dx[2]*(dx[3]^2 - 4*x[1]*d2x[1] + 3*x[2]*d2x[2] + x[3]*d2x[3]) - 8*x[1]*x[2]*d3x[1] - 4*x[1]^2*d3x[2] + 3*x[2]^2*d3x[2] + x[3]^2*d3x[2] + 2*x[2]*x[3]*d3x[3])/5.

@inline D3_M122(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (6*dx[1]^3 + 6*x[3]*dx[3]*d2x[1] - 24*x[1]*dx[2]*d2x[2] + 6*x[1]*dx[3]*d2x[3] + 6*dx[1]*(-4*dx[2]^2 + dx[3]^2 + 3*x[1]*d2x[1] - 4*x[2]*d2x[2] + x[3]*d2x[3]) + 3*x[1]^2*d3x[1] - 4*x[2]^2*d3x[1] + x[3]^2*d3x[1] - 8*x[2]*(3*dx[2]*d2x[1] + x[1]*d3x[2]) + 2*x[1]*x[3]*d3x[3])/5.

@inline D3_M113(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (-24*dx[1]^2*dx[3] + 6*dx[2]^2*dx[3] + 6*dx[3]^3 - 24*x[1]*dx[3]*d2x[1] + 6*x[2]*dx[3]*d2x[2] + 18*x[3]*dx[3]*d2x[3] - 24*dx[1]*(x[3]*d2x[1] + x[1]*d2x[3]) + 6*dx[2]*(x[3]*d2x[2] + x[2]*d2x[3]) - 8*x[1]*x[3]*d3x[1] + 2*x[2]*x[3]*d3x[2] + (-4*x[1]^2 + x[2]^2 + 3*x[3]^2)*d3x[3])/5.

@inline D3_M133(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (6*dx[1]^3 - 24*x[3]*dx[3]*d2x[1] + 6*x[1]*dx[2]*d2x[2] - 24*x[1]*dx[3]*d2x[3] + 6*dx[1]*(dx[2]^2 - 4*dx[3]^2 + 3*x[1]*d2x[1] + x[2]*d2x[2] - 4*x[3]*d2x[3]) + 3*x[1]^2*d3x[1] + x[2]^2*d3x[1] - 4*x[3]^2*d3x[1] + 2*x[2]*(3*dx[2]*d2x[1] + x[1]*d3x[2]) - 8*x[1]*x[3]*d3x[3])/5.

@inline D3_M123(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = -3*(x[3]*dx[2] + x[2]*dx[3])*d2x[1] - 3*dx[1]*(2*dx[2]*dx[3] + x[3]*d2x[2] + x[2]*d2x[3]) - x[2]*x[3]*d3x[1] - x[1]*(3*dx[3]*d2x[2] + 3*dx[2]*d2x[3] + x[3]*d3x[2] + x[2]*d3x[3])

@inline D3_M222(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (3*(6*dx[1]^2*dx[2] - 4*dx[2]^3 + 6*x[3]*dx[3]*d2x[2] + 6*dx[1]*(x[2]*d2x[1] + x[1]*d2x[2]) + 6*x[2]*dx[3]*d2x[3] + 6*dx[2]*(dx[3]^2 + x[1]*d2x[1] - 2*x[2]*d2x[2] + x[3]*d2x[3]) + 2*x[1]*x[2]*d3x[1] + x[1]^2*d3x[2] - 2*x[2]^2*d3x[2] + x[3]^2*d3x[2] + 2*x[2]*x[3]*d3x[3]))/5.

@inline D3_M223(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (6*dx[1]^2*dx[3] - 24*dx[2]^2*dx[3] + 6*dx[3]^3 + 6*x[1]*dx[3]*d2x[1] - 24*x[2]*dx[3]*d2x[2] + 18*x[3]*dx[3]*d2x[3] + 6*dx[1]*(x[3]*d2x[1] + x[1]*d2x[3]) - 24*dx[2]*(x[3]*d2x[2] + x[2]*d2x[3]) + 2*x[1]*x[3]*d3x[1] - 8*x[2]*x[3]*d3x[2] + (x[1]^2 - 4*x[2]^2 + 3*x[3]^2)*d3x[3])/5.

@inline D3_M233(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (6*dx[1]^2*dx[2] + 6*dx[2]^3 - 24*x[3]*dx[3]*d2x[2] + 6*dx[1]*(x[2]*d2x[1] + x[1]*d2x[2]) - 24*x[2]*dx[3]*d2x[3] + 6*dx[2]*(-4*dx[3]^2 + x[1]*d2x[1] + 3*x[2]*d2x[2] - 4*x[3]*d2x[3]) + 2*x[1]*x[2]*d3x[1] + x[1]^2*d3x[2] + 3*x[2]^2*d3x[2] - 4*x[3]^2*d3x[2] - 8*x[2]*x[3]*d3x[3])/5.

@inline D3_M333(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (6*(3*dx[1]^2*dx[3] + 3*dx[2]^2*dx[3] - 2*dx[3]^3 + 3*x[1]*dx[3]*d2x[1] + 3*x[2]*dx[3]*d2x[2] - 6*x[3]*dx[3]*d2x[3] + 3*dx[1]*(x[3]*d2x[1] + x[1]*d2x[3]) + 3*dx[2]*(x[3]*d2x[2] + x[2]*d2x[3]) + x[1]*x[3]*d3x[1] + x[2]*x[3]*d3x[2]) + 3*(x[1]^2 + x[2]^2 - 2*x[3]^2)*d3x[3])/5.

@inline D4_M1111(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (4*(-144*(dx[1]^2 + x[1]*d2x[1])*(dx[2]^2 + dx[3]^2 + x[2]*d2x[2] + x[3]*d2x[3]) + 18*(dx[2]^2 + dx[3]^2 + x[2]*d2x[2] + x[3]*d2x[3])^2 - 96*(x[2]*dx[2] + x[3]*dx[3])*(3*dx[1]*d2x[1] + x[1]*d3x[1]) - 96*x[1]*dx[1]*(3*dx[2]*d2x[2] + 3*dx[3]*d2x[3] + x[2]*d3x[2] + x[3]*d3x[3]) + 24*(x[2]*dx[2] + x[3]*dx[3])*(3*dx[2]*d2x[2] + 3*dx[3]*d2x[3] + x[2]*d3x[2] + x[3]*d3x[3]) - 12*(x[2]^2 + x[3]^2)*(3*d2x[1]^2 + 4*dx[1]*d3x[1] + x[1]*d4x[1]) + 8*(6*dx[1]^4 + 36*x[1]*dx[1]^2*d2x[1] + 12*x[1]^2*dx[1]*d3x[1] + x[1]^2*(9*d2x[1]^2 + x[1]*d4x[1])) - 12*x[1]^2*(3*d2x[2]^2 + 3*d2x[3]^2 + 4*dx[2]*d3x[2] + 4*dx[3]*d3x[3] + x[2]*d4x[2] + x[3]*d4x[3]) + 3*(x[2]^2 + x[3]^2)*(3*d2x[2]^2 + 3*d2x[3]^2 + 4*dx[2]*d3x[2] + 4*dx[3]*d3x[3] + x[2]*d4x[2] + x[3]*d4x[3])))/35.

@inline D4_M1112(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (96*dx[1]^3*dx[2] + 144*dx[1]^2*(x[2]*d2x[1] + x[1]*d2x[2]) - 12*dx[1]*(6*dx[2]^3 + 6*dx[2]*(dx[3]^2 - 4*x[1]*d2x[1] + 3*x[2]*d2x[2] + x[3]*d2x[3]) + 2*x[2]*(3*dx[3]*d2x[3] - 4*x[1]*d3x[1]) + (-4*x[1]^2 + 3*x[2]^2)*d3x[2] + x[3]^2*d3x[2] + 2*x[3]*(3*dx[3]*d2x[2] + x[2]*d3x[3])) - 6*(x[3]^2*(3*d2x[1]*d2x[2] + 2*dx[2]*d3x[1]) + 2*x[1]*(9*dx[2]^2*d2x[2] + 3*(dx[3]^2 - 2*x[1]*d2x[1])*d2x[2] + dx[2]*(6*dx[3]*d2x[3] - 4*x[1]*d3x[1])) + 2*x[3]*(x[1]*(3*d2x[2]*d2x[3] + 2*dx[3]*d3x[2]) + 2*dx[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3]))) - 3*x[2]^3*d4x[1] + x[1]*(4*x[1]^2 - 3*x[3]^2)*d4x[2] - 9*x[2]^2*(6*d2x[1]*d2x[2] + 4*dx[2]*d3x[1] + x[1]*d4x[2]) - 3*x[2]*(36*dx[2]^2*d2x[1] + 12*dx[3]^2*d2x[1] - 24*x[1]*d2x[1]^2 + 18*x[1]*d2x[2]^2 + 12*x[3]*d2x[1]*d2x[3] + 6*x[1]*d2x[3]^2 + 24*x[1]*dx[2]*d3x[2] + 8*dx[3]*(x[3]*d3x[1] + x[1]*d3x[3]) - 4*x[1]^2*d4x[1] + x[3]^2*d4x[1] + 2*x[1]*x[3]*d4x[3]))/7.

@inline D4_M1122(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (2*(12*dx[3]^4 + 324*(dx[1]^2 + x[1]*d2x[1])*(dx[2]^2 + x[2]*d2x[2]) + 72*x[3]*dx[3]^2*d2x[3] + 18*x[3]^2*d2x[3]^2 - 36*(dx[1]^2 + dx[2]^2 + x[1]*d2x[1] + x[2]*d2x[2])*(dx[3]^2 + x[3]*d2x[3]) + 216*x[2]*dx[2]*(3*dx[1]*d2x[1] + x[1]*d3x[1]) + 216*x[1]*dx[1]*(3*dx[2]*d2x[2] + x[2]*d3x[2]) - 24*x[3]*dx[3]*(3*dx[1]*d2x[1] + 3*dx[2]*d2x[2] + x[1]*d3x[1] + x[2]*d3x[2]) + 24*x[3]^2*dx[3]*d3x[3] - 24*(x[1]*dx[1] + x[2]*dx[2])*(3*dx[3]*d2x[3] + x[3]*d3x[3]) + 27*x[2]^2*(3*d2x[1]^2 + 4*dx[1]*d3x[1] + x[1]*d4x[1]) - 8*(6*dx[1]^4 + 36*x[1]*dx[1]^2*d2x[1] + 12*x[1]^2*dx[1]*d3x[1] + x[1]^2*(9*d2x[1]^2 + x[1]*d4x[1])) + 27*x[1]^2*(3*d2x[2]^2 + 4*dx[2]*d3x[2] + x[2]*d4x[2]) - 3*x[3]^2*(3*d2x[1]^2 + 3*d2x[2]^2 + 4*dx[1]*d3x[1] + 4*dx[2]*d3x[2] + x[1]*d4x[1] + x[2]*d4x[2]) - 8*(6*dx[2]^4 + 36*x[2]*dx[2]^2*d2x[2] + 12*x[2]^2*dx[2]*d3x[2] + x[2]^2*(9*d2x[2]^2 + x[2]*d4x[2])) + 2*x[3]^3*d4x[3] - 3*(x[1]^2 + x[2]^2)*(3*d2x[3]^2 + 4*dx[3]*d3x[3] + x[3]*d4x[3])))/35.

@inline D4_M1222(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-72*dx[1]^3*dx[2] - 108*dx[1]^2*(x[2]*d2x[1] + x[1]*d2x[2]) - 12*dx[1]*(-8*dx[2]^3 + 6*dx[2]*(dx[3]^2 + 3*x[1]*d2x[1] - 4*x[2]*d2x[2] + x[3]*d2x[3]) + 6*x[2]*(dx[3]*d2x[3] + x[1]*d3x[1]) + (3*x[1]^2 - 4*x[2]^2)*d3x[2] + x[3]^2*d3x[2] + 2*x[3]*(3*dx[3]*d2x[2] + x[2]*d3x[3])) - 6*(x[3]^2*(3*d2x[1]*d2x[2] + 2*dx[2]*d3x[1]) + 3*x[1]*(-8*dx[2]^2*d2x[2] + (2*dx[3]^2 + 3*x[1]*d2x[1])*d2x[2] + 2*dx[2]*(2*dx[3]*d2x[3] + x[1]*d3x[1])) + 2*x[3]*(x[1]*(3*d2x[2]*d2x[3] + 2*dx[3]*d3x[2]) + 2*dx[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3]))) + 4*x[2]^3*d4x[1] - 3*x[1]*(x[1]^2 + x[3]^2)*d4x[2] + 12*x[2]^2*(6*d2x[1]*d2x[2] + 4*dx[2]*d3x[1] + x[1]*d4x[2]) - 3*x[2]*(-48*dx[2]^2*d2x[1] + 12*dx[3]^2*d2x[1] + 18*x[1]*d2x[1]^2 - 24*x[1]*d2x[2]^2 + 12*x[3]*d2x[1]*d2x[3] + 6*x[1]*d2x[3]^2 - 32*x[1]*dx[2]*d3x[2] + 8*dx[3]*(x[3]*d3x[1] + x[1]*d3x[3]) + 3*x[1]^2*d4x[1] + x[3]^2*d4x[1] + 2*x[1]*x[3]*d4x[3]))/7.

@inline D4_M1113(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (96*dx[1]^3*dx[3] + 144*dx[1]^2*(x[3]*d2x[1] + x[1]*d2x[3]) - 12*dx[1]*(6*dx[2]^2*dx[3] + 6*dx[3]^3 + 6*dx[2]*(x[3]*d2x[2] + x[2]*d2x[3]) + 6*dx[3]*(-4*x[1]*d2x[1] + x[2]*d2x[2] + 3*x[3]*d2x[3]) + 2*x[3]*(-4*x[1]*d3x[1] + x[2]*d3x[2]) + (-4*x[1]^2 + x[2]^2 + 3*x[3]^2)*d3x[3]) - 6*(x[2]^2*(3*d2x[1]*d2x[3] + 2*dx[3]*d3x[1]) + 2*x[1]*(6*dx[2]*dx[3]*d2x[2] + 3*dx[2]^2*d2x[3] + 3*(3*dx[3]^2 - 2*x[1]*d2x[1])*d2x[3] - 4*x[1]*dx[3]*d3x[1]) + 2*x[2]*(x[1]*(3*d2x[2]*d2x[3] + 2*dx[3]*d3x[2]) + 2*dx[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3]))) - 3*x[3]^3*d4x[1] - 3*x[3]*(12*dx[2]^2*d2x[1] + 36*dx[3]^2*d2x[1] - 24*x[1]*d2x[1]^2 + 12*x[2]*d2x[1]*d2x[2] + 6*x[1]*d2x[2]^2 + 18*x[1]*d2x[3]^2 + 8*dx[2]*(x[2]*d3x[1] + x[1]*d3x[2]) + 24*x[1]*dx[3]*d3x[3] - 4*x[1]^2*d4x[1] + x[2]^2*d4x[1] + 2*x[1]*x[2]*d4x[2]) + x[1]*(4*x[1]^2 - 3*x[2]^2)*d4x[3] - 9*x[3]^2*(6*d2x[1]*d2x[3] + 4*dx[3]*d3x[1] + x[1]*d4x[3]))/7.

@inline D4_M1133(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (2*(12*dx[2]^4 + 72*x[2]*dx[2]^2*d2x[2] + 18*x[2]^2*d2x[2]^2 - 36*(dx[1]^2 + x[1]*d2x[1])*(dx[2]^2 - 9*dx[3]^2 + x[2]*d2x[2] - 9*x[3]*d2x[3]) - 36*(dx[2]^2 + x[2]*d2x[2])*(dx[3]^2 + x[3]*d2x[3]) - 24*(x[2]*dx[2] - 9*x[3]*dx[3])*(3*dx[1]*d2x[1] + x[1]*d3x[1]) + 24*x[2]^2*dx[2]*d3x[2] - 24*x[3]*dx[3]*(3*dx[2]*d2x[2] + x[2]*d3x[2]) - 24*x[1]*dx[1]*(3*dx[2]*d2x[2] - 27*dx[3]*d2x[3] + x[2]*d3x[2] - 9*x[3]*d3x[3]) - 24*x[2]*dx[2]*(3*dx[3]*d2x[3] + x[3]*d3x[3]) - 3*(x[2]^2 - 9*x[3]^2)*(3*d2x[1]^2 + 4*dx[1]*d3x[1] + x[1]*d4x[1]) - 8*(6*dx[1]^4 + 36*x[1]*dx[1]^2*d2x[1] + 12*x[1]^2*dx[1]*d3x[1] + x[1]^2*(9*d2x[1]^2 + x[1]*d4x[1])) + 2*x[2]^3*d4x[2] - 3*x[3]^2*(3*d2x[2]^2 + 4*dx[2]*d3x[2] + x[2]*d4x[2]) - 3*x[1]^2*(3*d2x[2]^2 - 27*d2x[3]^2 + 4*dx[2]*d3x[2] - 36*dx[3]*d3x[3] + x[2]*d4x[2] - 9*x[3]*d4x[3]) - 3*x[2]^2*(3*d2x[3]^2 + 4*dx[3]*d3x[3] + x[3]*d4x[3]) - 8*(6*dx[3]^4 + 36*x[3]*dx[3]^2*d2x[3] + 12*x[3]^2*dx[3]*d3x[3] + x[3]^2*(9*d2x[3]^2 + x[3]*d4x[3]))))/35.

@inline D4_M1333(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-72*dx[1]^3*dx[3] - 108*dx[1]^2*(x[3]*d2x[1] + x[1]*d2x[3]) + 12*dx[1]*(-6*dx[2]^2*dx[3] + 8*dx[3]^3 - 6*dx[2]*(x[3]*d2x[2] + x[2]*d2x[3]) - 6*dx[3]*(3*x[1]*d2x[1] + x[2]*d2x[2] - 4*x[3]*d2x[3]) - 2*x[3]*(3*x[1]*d3x[1] + x[2]*d3x[2]) - (3*x[1]^2 + x[2]^2 - 4*x[3]^2)*d3x[3]) - 6*(x[2]^2*(3*d2x[1]*d2x[3] + 2*dx[3]*d3x[1]) + 3*x[1]*(4*dx[2]*dx[3]*d2x[2] + 2*dx[2]^2*d2x[3] - 8*dx[3]^2*d2x[3] + 3*x[1]*d2x[1]*d2x[3] + 2*x[1]*dx[3]*d3x[1]) + 2*x[2]*(x[1]*(3*d2x[2]*d2x[3] + 2*dx[3]*d3x[2]) + 2*dx[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3]))) + 4*x[3]^3*d4x[1] - 3*x[3]*(12*dx[2]^2*d2x[1] - 48*dx[3]^2*d2x[1] + 18*x[1]*d2x[1]^2 + 12*x[2]*d2x[1]*d2x[2] + 6*x[1]*d2x[2]^2 - 24*x[1]*d2x[3]^2 + 8*dx[2]*(x[2]*d3x[1] + x[1]*d3x[2]) - 32*x[1]*dx[3]*d3x[3] + 3*x[1]^2*d4x[1] + x[2]^2*d4x[1] + 2*x[1]*x[2]*d4x[2]) - 3*x[1]*(x[1]^2 + x[2]^2)*d4x[3] + 12*x[3]^2*(6*d2x[1]*d2x[3] + 4*dx[3]*d3x[1] + x[1]*d4x[3]))/7.

@inline D4_M1123(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-24*dx[2]^3*dx[3] + 36*x[2]*x[3]*d2x[1]^2 - 36*x[3]*dx[3]^2*d2x[2] + 72*x[1]*x[3]*d2x[1]*d2x[2] - 18*x[2]*x[3]*d2x[2]^2 - 36*x[2]*dx[3]^2*d2x[3] + 72*x[1]*x[2]*d2x[1]*d2x[3] + 36*x[1]^2*d2x[2]*d2x[3] - 18*x[2]^2*d2x[2]*d2x[3] - 18*x[3]^2*d2x[2]*d2x[3] - 18*x[2]*x[3]*d2x[3]^2 - 36*dx[2]^2*(x[3]*d2x[2] + x[2]*d2x[3]) + 72*dx[1]^2*(2*dx[2]*dx[3] + x[3]*d2x[2] + x[2]*d2x[3]) + 48*x[1]*x[2]*dx[3]*d3x[1] + 24*x[1]^2*dx[3]*d3x[2] - 12*x[2]^2*dx[3]*d3x[2] - 12*x[3]^2*dx[3]*d3x[2] - 24*x[2]*x[3]*dx[3]*d3x[3] - 12*dx[2]*(2*(dx[3]^3 + 3*dx[3]*(-2*x[1]*d2x[1] + x[2]*d2x[2] + x[3]*d2x[3]) + x[3]*(-2*x[1]*d3x[1] + x[2]*d3x[2])) + (-2*x[1]^2 + x[2]^2 + x[3]^2)*d3x[3]) + 48*dx[1]*(3*x[1]*(dx[3]*d2x[2] + dx[2]*d2x[3]) + x[3]*(3*dx[2]*d2x[1] + x[2]*d3x[1] + x[1]*d3x[2]) + x[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3])) + 12*x[1]*x[2]*x[3]*d4x[1] + 6*x[1]^2*x[3]*d4x[2] - 3*x[2]^2*x[3]*d4x[2] - x[3]^3*d4x[2] - x[2]*(-6*x[1]^2 + x[2]^2 + 3*x[3]^2)*d4x[3])/7.

@inline D4_M1223(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-24*dx[1]^3*dx[3] - 36*dx[1]^2*(x[3]*d2x[1] + x[1]*d2x[3]) + 12*dx[1]*(12*dx[2]^2*dx[3] + 12*dx[2]*(x[3]*d2x[2] + x[2]*d2x[3]) - 2*(dx[3]^3 + 3*dx[3]*(x[1]*d2x[1] - 2*x[2]*d2x[2] + x[3]*d2x[3]) + x[3]*(x[1]*d3x[1] - 2*x[2]*d3x[2])) - (x[1]^2 - 2*x[2]^2 + x[3]^2)*d3x[3]) + 6*(x[2]^2*(6*d2x[1]*d2x[3] + 4*dx[3]*d3x[1]) + x[1]*(24*dx[2]*dx[3]*d2x[2] + 12*dx[2]^2*d2x[3] - 3*(2*dx[3]^2 + x[1]*d2x[1])*d2x[3] - 2*x[1]*dx[3]*d3x[1]) + 4*x[2]*(x[1]*(3*d2x[2]*d2x[3] + 2*dx[3]*d3x[2]) + 2*dx[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3]))) - x[3]^3*d4x[1] + 3*x[3]*(24*dx[2]^2*d2x[1] - 12*dx[3]^2*d2x[1] - 6*x[1]*d2x[1]^2 + 24*x[2]*d2x[1]*d2x[2] + 12*x[1]*d2x[2]^2 - 6*x[1]*d2x[3]^2 + 16*dx[2]*(x[2]*d3x[1] + x[1]*d3x[2]) - 8*x[1]*dx[3]*d3x[3] - x[1]^2*d4x[1] + 2*x[2]^2*d4x[1] + 4*x[1]*x[2]*d4x[2]) - x[1]*(x[1]^2 - 6*x[2]^2)*d4x[3] - 3*x[3]^2*(6*d2x[1]*d2x[3] + 4*dx[3]*d3x[1] + x[1]*d4x[3]))/7.

@inline D4_M1233(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-24*dx[1]^3*dx[2] - 36*dx[1]^2*(x[2]*d2x[1] + x[1]*d2x[2]) - 12*dx[1]*(2*dx[2]^3 + 6*dx[2]*(-2*dx[3]^2 + x[1]*d2x[1] + x[2]*d2x[2] - 2*x[3]*d2x[3]) + 2*x[2]*(-6*dx[3]*d2x[3] + x[1]*d3x[1]) + (x[1]^2 + x[2]^2)*d3x[2] - 2*x[3]^2*d3x[2] - 4*x[3]*(3*dx[3]*d2x[2] + x[2]*d3x[3])) + 6*(x[3]^2*(6*d2x[1]*d2x[2] + 4*dx[2]*d3x[1]) + x[1]*(-6*dx[2]^2*d2x[2] + 3*(4*dx[3]^2 - x[1]*d2x[1])*d2x[2] + dx[2]*(24*dx[3]*d2x[3] - 2*x[1]*d3x[1])) + 4*x[3]*(x[1]*(3*d2x[2]*d2x[3] + 2*dx[3]*d3x[2]) + 2*dx[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3]))) - x[2]^3*d4x[1] - x[1]*(x[1]^2 - 6*x[3]^2)*d4x[2] - 3*x[2]^2*(6*d2x[1]*d2x[2] + 4*dx[2]*d3x[1] + x[1]*d4x[2]) + 3*x[2]*(-12*dx[2]^2*d2x[1] + 24*dx[3]^2*d2x[1] - 6*x[1]*d2x[1]^2 - 6*x[1]*d2x[2]^2 + 24*x[3]*d2x[1]*d2x[3] + 12*x[1]*d2x[3]^2 - 8*x[1]*dx[2]*d3x[2] + 16*dx[3]*(x[3]*d3x[1] + x[1]*d3x[3]) - x[1]^2*d4x[1] + 2*x[3]^2*d4x[1] + 4*x[1]*x[3]*d4x[3]))/7.

@inline D4_M2222(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (4*(210*dx[2]^4 + 1260*x[2]*dx[2]^2*d2x[2] + 315*x[2]^2*d2x[2]^2 - 180*(dx[2]^2 + x[2]*d2x[2])*(dx[1]^2 + dx[2]^2 + dx[3]^2 + x[1]*d2x[1] + x[2]*d2x[2] + x[3]*d2x[3]) + 420*x[2]^2*dx[2]*d3x[2] - 120*(x[1]*dx[1] + x[2]*dx[2] + x[3]*dx[3])*(3*dx[2]*d2x[2] + x[2]*d3x[2]) - 120*x[2]*dx[2]*(3*dx[1]*d2x[1] + 3*dx[2]*d2x[2] + 3*dx[3]*d2x[3] + x[1]*d3x[1] + x[2]*d3x[2] + x[3]*d3x[3]) + 35*x[2]^3*d4x[2] - 15*(x[1]^2 + x[2]^2 + x[3]^2)*(3*d2x[2]^2 + 4*dx[2]*d3x[2] + x[2]*d4x[2]) - 15*x[2]^2*(3*d2x[1]^2 + 3*d2x[2]^2 + 3*d2x[3]^2 + 4*dx[1]*d3x[1] + 4*dx[2]*d3x[2] + 4*dx[3]*d3x[3] + x[1]*d4x[1] + x[2]*d4x[2] + x[3]*d4x[3]) + 3*(6*(dx[1]^2 + dx[2]^2 + dx[3]^2 + x[1]*d2x[1] + x[2]*d2x[2] + x[3]*d2x[3])^2 + 8*(x[1]*dx[1] + x[2]*dx[2] + x[3]*dx[3])*(3*dx[1]*d2x[1] + 3*dx[2]*d2x[2] + 3*dx[3]*d2x[3] + x[1]*d3x[1] + x[2]*d3x[2] + x[3]*d3x[3]) + (x[1]^2 + x[2]^2 + x[3]^2)*(3*d2x[1]^2 + 3*d2x[2]^2 + 3*d2x[3]^2 + 4*dx[1]*d3x[1] + 4*dx[2]*d3x[2] + 4*dx[3]*d3x[3] + x[1]*d4x[1] + x[2]*d4x[2] + x[3]*d4x[3]))))/35.

@inline D4_M2223(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (96*dx[2]^3*dx[3] + 144*dx[2]^2*(x[3]*d2x[2] + x[2]*d2x[3]) - 36*dx[1]^2*(2*dx[2]*dx[3] + x[3]*d2x[2] + x[2]*d2x[3]) + 12*dx[2]*(-6*dx[3]^3 - 6*dx[3]*(x[1]*d2x[1] - 4*x[2]*d2x[2] + 3*x[3]*d2x[3]) - 2*x[3]*(x[1]*d3x[1] - 4*x[2]*d3x[2]) - (x[1]^2 - 4*x[2]^2 + 3*x[3]^2)*d3x[3]) - 24*dx[1]*(3*x[1]*(dx[3]*d2x[2] + dx[2]*d2x[3]) + x[3]*(3*dx[2]*d2x[1] + x[2]*d3x[1] + x[1]*d3x[2]) + x[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3])) - 3*(6*d2x[2]*(2*x[3]*(3*dx[3]^2 + x[1]*d2x[1]) + (x[1]^2 + 3*x[3]^2)*d2x[3]) + 4*(x[1]^2 + 3*x[3]^2)*dx[3]*d3x[2] + 2*x[2]*(6*(3*dx[3]^2 + x[1]*d2x[1])*d2x[3] + 4*x[1]*dx[3]*d3x[1] + x[3]*(3*(d2x[1]^2 - 4*d2x[2]^2 + 3*d2x[3]^2 + 4*dx[3]*d3x[3]) + x[1]*d4x[1])) + x[3]*(x[1]^2 + x[3]^2)*d4x[2] - 4*x[2]^2*(6*d2x[2]*d2x[3] + 4*dx[3]*d3x[2] + x[3]*d4x[2])) + x[2]*(-3*x[1]^2 + 4*x[2]^2 - 9*x[3]^2)*d4x[3])/7.

@inline D4_M2233(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (2*(12*dx[1]^4 + 72*x[1]*dx[1]^2*d2x[1] + 18*x[1]^2*d2x[1]^2 + 324*(dx[2]^2 + x[2]*d2x[2])*(dx[3]^2 + x[3]*d2x[3]) - 36*(dx[1]^2 + x[1]*d2x[1])*(dx[2]^2 + dx[3]^2 + x[2]*d2x[2] + x[3]*d2x[3]) + 24*x[1]^2*dx[1]*d3x[1] - 24*(x[2]*dx[2] + x[3]*dx[3])*(3*dx[1]*d2x[1] + x[1]*d3x[1]) + 216*x[3]*dx[3]*(3*dx[2]*d2x[2] + x[2]*d3x[2]) + 216*x[2]*dx[2]*(3*dx[3]*d2x[3] + x[3]*d3x[3]) - 24*x[1]*dx[1]*(3*dx[2]*d2x[2] + 3*dx[3]*d2x[3] + x[2]*d3x[2] + x[3]*d3x[3]) + 2*x[1]^3*d4x[1] - 3*(x[2]^2 + x[3]^2)*(3*d2x[1]^2 + 4*dx[1]*d3x[1] + x[1]*d4x[1]) + 27*x[3]^2*(3*d2x[2]^2 + 4*dx[2]*d3x[2] + x[2]*d4x[2]) - 8*(6*dx[2]^4 + 36*x[2]*dx[2]^2*d2x[2] + 12*x[2]^2*dx[2]*d3x[2] + x[2]^2*(9*d2x[2]^2 + x[2]*d4x[2])) + 27*x[2]^2*(3*d2x[3]^2 + 4*dx[3]*d3x[3] + x[3]*d4x[3]) - 3*x[1]^2*(3*d2x[2]^2 + 3*d2x[3]^2 + 4*dx[2]*d3x[2] + 4*dx[3]*d3x[3] + x[2]*d4x[2] + x[3]*d4x[3]) - 8*(6*dx[3]^4 + 36*x[3]*dx[3]^2*d2x[3] + 12*x[3]^2*dx[3]*d3x[3] + x[3]^2*(9*d2x[3]^2 + x[3]*d4x[3]))))/35.

@inline D4_M2333(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-72*dx[2]^3*dx[3] - 18*x[2]*x[3]*d2x[1]^2 + 144*x[3]*dx[3]^2*d2x[2] - 36*x[1]*x[3]*d2x[1]*d2x[2] - 54*x[2]*x[3]*d2x[2]^2 + 144*x[2]*dx[3]^2*d2x[3] - 36*x[1]*x[2]*d2x[1]*d2x[3] - 18*x[1]^2*d2x[2]*d2x[3] - 54*x[2]^2*d2x[2]*d2x[3] + 72*x[3]^2*d2x[2]*d2x[3] + 72*x[2]*x[3]*d2x[3]^2 - 108*dx[2]^2*(x[3]*d2x[2] + x[2]*d2x[3]) - 36*dx[1]^2*(2*dx[2]*dx[3] + x[3]*d2x[2] + x[2]*d2x[3]) - 24*x[1]*x[2]*dx[3]*d3x[1] - 12*x[1]^2*dx[3]*d3x[2] - 36*x[2]^2*dx[3]*d3x[2] + 48*x[3]^2*dx[3]*d3x[2] + 96*x[2]*x[3]*dx[3]*d3x[3] - 12*dx[2]*(-8*dx[3]^3 + 6*dx[3]*(x[1]*d2x[1] + 3*x[2]*d2x[2] - 4*x[3]*d2x[3]) + 2*x[3]*(x[1]*d3x[1] + 3*x[2]*d3x[2]) + (x[1]^2 + 3*x[2]^2 - 4*x[3]^2)*d3x[3]) - 24*dx[1]*(3*x[1]*(dx[3]*d2x[2] + dx[2]*d2x[3]) + x[3]*(3*dx[2]*d2x[1] + x[2]*d3x[1] + x[1]*d3x[2]) + x[2]*(3*dx[3]*d2x[1] + x[1]*d3x[3])) - 6*x[1]*x[2]*x[3]*d4x[1] - 3*x[1]^2*x[3]*d4x[2] - 9*x[2]^2*x[3]*d4x[2] + 4*x[3]^3*d4x[2] - 3*x[2]*(x[1]^2 + x[2]^2 - 4*x[3]^2)*d4x[3])/7.

@inline D4_M3333(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (4*(18*(dx[1]^2 + dx[2]^2 + x[1]*d2x[1] + x[2]*d2x[2])^2 - 144*(dx[1]^2 + dx[2]^2 + x[1]*d2x[1] + x[2]*d2x[2])*(dx[3]^2 + x[3]*d2x[3]) + 24*(x[1]*dx[1] + x[2]*dx[2])*(3*dx[1]*d2x[1] + 3*dx[2]*d2x[2] + x[1]*d3x[1] + x[2]*d3x[2]) - 96*x[3]*dx[3]*(3*dx[1]*d2x[1] + 3*dx[2]*d2x[2] + x[1]*d3x[1] + x[2]*d3x[2]) - 96*(x[1]*dx[1] + x[2]*dx[2])*(3*dx[3]*d2x[3] + x[3]*d3x[3]) + 3*(x[1]^2 + x[2]^2)*(3*d2x[1]^2 + 3*d2x[2]^2 + 4*dx[1]*d3x[1] + 4*dx[2]*d3x[2] + x[1]*d4x[1] + x[2]*d4x[2]) - 12*x[3]^2*(3*d2x[1]^2 + 3*d2x[2]^2 + 4*dx[1]*d3x[1] + 4*dx[2]*d3x[2] + x[1]*d4x[1] + x[2]*d4x[2]) - 12*(x[1]^2 + x[2]^2)*(3*d2x[3]^2 + 4*dx[3]*d3x[3] + x[3]*d4x[3]) + 8*(6*dx[3]^4 + 36*x[3]*dx[3]^2*d2x[3] + 12*x[3]^2*dx[3]*d3x[3] + x[3]^2*(9*d2x[3]^2 + x[3]*d4x[3]))))/35.

@inline D2_S12(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (3*(x[1]*dx[1] - x[2]*dx[2])*d2x[3] + x[3]*(-3*dx[1]*d2x[1] + 3*dx[2]*d2x[2] - x[1]*d3x[1] + x[2]*d3x[2]) + (x[1]^2 - x[2]^2)*d3x[3])/2.

@inline D2_S13(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (3*(-(x[1]*dx[1]) + x[3]*dx[3])*d2x[2] + (-x[1]^2 + x[3]^2)*d3x[2] + x[2]*(3*dx[1]*d2x[1] - 3*dx[3]*d2x[3] + x[1]*d3x[1] - x[3]*d3x[3]))/2.

@inline D2_S23(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (3*x[1]*(-(dx[2]*d2x[2]) + dx[3]*d2x[3]) + x[2]^2*d3x[1] - x[3]^2*d3x[1] + x[2]*(3*dx[2]*d2x[1] - x[1]*d3x[2]) + x[3]*(-3*dx[3]*d2x[1] + x[1]*d3x[3]))/2.

@inline D2_S11(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (x[3]*dx[2] - x[2]*dx[3])*d2x[1] + 2*dx[1]*(x[3]*d2x[2] - x[2]*d2x[3]) + x[1]*(dx[3]*d2x[2] - dx[2]*d2x[3] + x[3]*d3x[2] - x[2]*d3x[3])

@inline D2_S22(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = (-(x[3]*dx[1]) + x[1]*dx[3])*d2x[2] + 2*dx[2]*(-(x[3]*d2x[1]) + x[1]*d2x[3]) + x[2]*(-(dx[3]*d2x[1]) + dx[1]*d2x[3] - x[3]*d3x[1] + x[1]*d3x[3])

@inline D2_S33(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64})::Float64 = 2*dx[3]*(x[2]*d2x[1] - x[1]*d2x[2]) + (x[2]*dx[1] - x[1]*dx[2])*d2x[3] + x[3]*(dx[2]*d2x[1] - dx[1]*d2x[2] + x[2]*d3x[1] - x[1]*d3x[2])

@inline D3_S111(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (6*(x[3]*d2x[2] - x[2]*d2x[3])*(-4*dx[1]^2 + dx[2]^2 + dx[3]^2 - 4*x[1]*d2x[1] + x[2]*d2x[2] + x[3]*d2x[3]) + 6*(4*x[1]*dx[1] - x[2]*dx[2] - x[3]*dx[3])*(-(dx[3]*d2x[2]) + dx[2]*d2x[3] - x[3]*d3x[2] + x[2]*d3x[3]) + 2*(x[3]*dx[2] - x[2]*dx[3])*(-12*dx[1]*d2x[1] + 3*dx[2]*d2x[2] + 3*dx[3]*d2x[3] - 4*x[1]*d3x[1] + x[2]*d3x[2] + x[3]*d3x[3]) + (4*x[1]^2 - x[2]^2 - x[3]^2)*(-2*dx[3]*d3x[2] + 2*dx[2]*d3x[3] - x[3]*d4x[2] + x[2]*d4x[3]))/5.

@inline D3_S112(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (4*x[2]*(x[1]*(9*d2x[2]*d2x[3] - 2*dx[3]*d3x[2]) + 2*dx[2]*(3*dx[3]*d2x[1] + 12*dx[1]*d2x[3] + 7*x[1]*d3x[3])) - x[3]^3*d4x[1] + x[3]*(48*dx[1]^2*d2x[1] - 36*dx[2]^2*d2x[1] - 12*dx[3]^2*d2x[1] + 24*x[1]*d2x[1]^2 - 36*x[2]*d2x[1]*d2x[2] - 30*x[1]*d2x[2]^2 + 6*x[1]*d2x[3]^2 - 8*dx[2]*(2*x[2]*d3x[1] + 5*x[1]*d3x[2]) - 32*dx[1]*(3*dx[2]*d2x[2] - x[1]*d3x[1] + x[2]*d3x[2]) + 8*x[1]*dx[3]*d3x[3] + 4*x[1]^2*d4x[1] - x[2]^2*d4x[1] - 10*x[1]*x[2]*d4x[2]) + x[3]^2*(-6*d2x[1]*d2x[3] - 8*dx[3]*d3x[1] + x[1]*d4x[3]) + x[2]^2*(30*d2x[1]*d2x[3] + 8*dx[3]*d3x[1] + 32*dx[1]*d3x[3] + 11*x[1]*d4x[3]) - 4*x[1]*(6*dx[2]*dx[3]*d2x[2] + 12*dx[1]^2*d2x[3] - 9*dx[2]^2*d2x[3] - 3*(dx[3]^2 - 2*x[1]*d2x[1])*d2x[3] + 8*x[1]*dx[1]*d3x[3] + x[1]^2*d4x[3]))/15.

@inline D3_S122(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (8*x[2]^2*(3*d2x[2]*d2x[3] + 4*dx[2]*d3x[3]) - 2*x[1]*(12*dx[1]*(dx[3]*d2x[2] + 4*dx[2]*d2x[3]) + x[1]*(15*d2x[2]*d2x[3] + 4*dx[3]*d3x[2] + 16*dx[2]*d3x[3])) + x[3]^3*d4x[2] + x[3]*(36*dx[1]^2*d2x[2] + 12*(-4*dx[2]^2 + dx[3]^2)*d2x[2] + 8*dx[1]*(12*dx[2]*d2x[1] + 5*x[2]*d3x[1] + 2*x[1]*d3x[2]) + 2*x[2]*(15*d2x[1]^2 - 12*d2x[2]^2 - 3*d2x[3]^2 - 16*dx[2]*d3x[2] - 4*dx[3]*d3x[3] + 5*x[1]*d4x[1]) - 4*x[2]^2*d4x[2] + x[1]*(36*d2x[1]*d2x[2] + 32*dx[2]*d3x[1] + x[1]*d4x[2])) + 4*x[2]^3*d4x[3] - x[2]*(4*(9*dx[1]^2*d2x[3] + 3*(-4*dx[2]^2 + dx[3]^2 + 3*x[1]*d2x[1])*d2x[3] - 2*x[1]*dx[3]*d3x[1] + dx[1]*(-6*dx[3]*d2x[1] + 14*x[1]*d3x[3])) + 11*x[1]^2*d4x[3]) + x[3]^2*(6*d2x[2]*d2x[3] + 8*dx[3]*d3x[2] - x[2]*d4x[3]))/15.

@inline D3_S113(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-2*x[3]^2*(15*d2x[1]*d2x[2] + 4*dx[2]*d3x[1] + 16*dx[1]*d3x[2]) + 4*x[1]*(12*dx[1]^2*d2x[2] - 3*(dx[2]^2 + 3*dx[3]^2 - 2*x[1]*d2x[1])*d2x[2] + 6*dx[2]*dx[3]*d2x[3] + 8*x[1]*dx[1]*d3x[2]) - 4*x[3]*(3*d2x[2]*(8*dx[1]*dx[3] + 3*x[1]*d2x[3]) + 14*x[1]*dx[3]*d3x[2] + dx[2]*(6*dx[3]*d2x[1] - 2*x[1]*d3x[3])) + x[2]^3*d4x[1] + x[1]*(4*x[1]^2 - 11*x[3]^2)*d4x[2] + x[2]^2*(6*d2x[1]*d2x[2] + 8*dx[2]*d3x[1] - x[1]*d4x[2]) + x[2]*(-48*dx[1]^2*d2x[1] + 12*dx[2]^2*d2x[1] + 36*dx[3]^2*d2x[1] - 24*x[1]*d2x[1]^2 - 6*x[1]*d2x[2]^2 + 36*x[3]*d2x[1]*d2x[3] + 30*x[1]*d2x[3]^2 + 16*x[3]*dx[3]*d3x[1] - 8*x[1]*dx[2]*d3x[2] + 40*x[1]*dx[3]*d3x[3] + 32*dx[1]*(3*dx[3]*d2x[3] - x[1]*d3x[1] + x[3]*d3x[3]) - 4*x[1]^2*d4x[1] + x[3]^2*d4x[1] + 10*x[1]*x[3]*d4x[3]))/15.

@inline D3_S133(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-2*x[2]^2*(3*d2x[2]*d2x[3] + 4*dx[2]*d3x[3]) - 4*x[2]*(9*dx[1]^2*d2x[3] + 3*(dx[2]^2 - 4*dx[3]^2 + 3*x[1]*d2x[1])*d2x[3] + 8*x[1]*dx[3]*d3x[1] + 4*dx[1]*(6*dx[3]*d2x[1] + x[1]*d3x[3])) + 2*x[1]*(12*dx[1]*(4*dx[3]*d2x[2] + dx[2]*d2x[3]) + x[1]*(15*d2x[2]*d2x[3] + 16*dx[3]*d3x[2] + 4*dx[2]*d3x[3])) - 4*x[3]^3*d4x[2] + x[3]*(36*dx[1]^2*d2x[2] + 12*(dx[2]^2 - 4*dx[3]^2)*d2x[2] - 8*dx[1]*(3*dx[2]*d2x[1] + 5*x[2]*d3x[1] - 7*x[1]*d3x[2]) + x[2]*(-30*d2x[1]^2 + 6*d2x[2]^2 + 24*d2x[3]^2 + 8*dx[2]*d3x[2] + 32*dx[3]*d3x[3] - 10*x[1]*d4x[1]) + x[2]^2*d4x[2] + x[1]*(36*d2x[1]*d2x[2] - 8*dx[2]*d3x[1] + 11*x[1]*d4x[2])) - x[2]*(x[1]^2 + x[2]^2)*d4x[3] + 4*x[3]^2*(-6*d2x[2]*d2x[3] - 8*dx[3]*d3x[2] + x[2]*d4x[3]))/15.

@inline D3_S123(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (x[1]*(12*dx[1]*(dx[2]*d2x[2] - dx[3]*d2x[3]) + x[1]*(3*d2x[2]^2 - 3*d2x[3]^2 + 4*dx[2]*d3x[2] - 4*dx[3]*d3x[3])) + x[3]^2*(3*d2x[1]^2 - 3*d2x[2]^2 + 4*dx[1]*d3x[1] - 4*dx[2]*d3x[2] + x[1]*d4x[1]) + x[2]*(-4*x[3]*dx[3]*d3x[2] + 4*dx[1]*(-3*dx[2]*d2x[1] + x[1]*d3x[2]) + 4*dx[2]*(3*dx[3]*d2x[3] - x[1]*d3x[1] + x[3]*d3x[3]) + (x[1]^2 - x[3]^2)*d4x[2]) + x[3]*(4*dx[3]*(-3*dx[2]*d2x[2] + x[1]*d3x[1]) + 4*dx[1]*(3*dx[3]*d2x[1] - x[1]*d3x[3]) - x[1]^2*d4x[3]) + x[2]^2*(-3*d2x[1]^2 + 3*d2x[3]^2 - 4*dx[1]*d3x[1] + 4*dx[3]*d3x[3] - x[1]*d4x[1] + x[3]*d4x[3]))/3.

@inline D3_S222(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (-6*(x[3]*d2x[1] - x[1]*d2x[3])*(dx[1]^2 - 4*dx[2]^2 + dx[3]^2 + x[1]*d2x[1] - 4*x[2]*d2x[2] + x[3]*d2x[3]) + 6*(x[1]*dx[1] - 4*x[2]*dx[2] + x[3]*dx[3])*(-(dx[3]*d2x[1]) + dx[1]*d2x[3] - x[3]*d3x[1] + x[1]*d3x[3]) - 2*(x[3]*dx[1] - x[1]*dx[3])*(3*dx[1]*d2x[1] - 12*dx[2]*d2x[2] + 3*dx[3]*d2x[3] + x[1]*d3x[1] - 4*x[2]*d3x[2] + x[3]*d3x[3]) + (x[1]^2 - 4*x[2]^2 + x[3]^2)*(-2*dx[3]*d3x[1] + 2*dx[1]*d3x[3] - x[3]*d4x[1] + x[1]*d4x[3]))/5.

@inline D3_S223(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (2*x[3]^2*(15*d2x[1]*d2x[2] + 16*dx[2]*d3x[1] + 4*dx[1]*d3x[2]) - 2*x[1]*(6*dx[1]^2*d2x[2] + 3*(-8*dx[2]^2 + 6*dx[3]^2 + x[1]*d2x[1])*d2x[2] + 48*dx[2]*dx[3]*d2x[3] + 4*x[1]*dx[1]*d3x[2]) - 4*x[3]*(-6*dx[1]*dx[3]*d2x[2] + 9*x[1]*d2x[2]*d2x[3] + 4*x[1]*dx[3]*d3x[2] + 8*dx[2]*(-3*dx[3]*d2x[1] + x[1]*d3x[3])) - 4*x[2]^3*d4x[1] - x[1]*(x[1]^2 + x[3]^2)*d4x[2] + 4*x[2]^2*(-6*d2x[1]*d2x[2] - 8*dx[2]*d3x[1] + x[1]*d4x[2]) + x[2]*(12*dx[1]^2*d2x[1] - 48*dx[2]^2*d2x[1] + 36*dx[3]^2*d2x[1] + 6*x[1]*d2x[1]^2 + 24*x[1]*d2x[2]^2 + 36*x[3]*d2x[1]*d2x[3] - 30*x[1]*d2x[3]^2 + 56*x[3]*dx[3]*d3x[1] + 32*x[1]*dx[2]*d3x[2] - 40*x[1]*dx[3]*d3x[3] - 8*dx[1]*(3*dx[3]*d2x[3] - x[1]*d3x[1] + x[3]*d3x[3]) + x[1]^2*d4x[1] + 11*x[3]^2*d4x[1] - 10*x[1]*x[3]*d4x[3]))/15.

@inline D3_S233(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (4*x[2]*(x[1]*(9*d2x[2]*d2x[3] + 8*dx[3]*d3x[2]) + dx[2]*(-24*dx[3]*d2x[1] - 6*dx[1]*d2x[3] + 4*x[1]*d3x[3])) + 4*x[3]^3*d4x[1] - x[3]*(12*dx[1]^2*d2x[1] + 36*dx[2]^2*d2x[1] - 48*dx[3]^2*d2x[1] + 6*x[1]*d2x[1]^2 + 36*x[2]*d2x[1]*d2x[2] - 30*x[1]*d2x[2]^2 + 24*x[1]*d2x[3]^2 + 8*dx[2]*(7*x[2]*d3x[1] - 5*x[1]*d3x[2]) - 8*dx[1]*(3*dx[2]*d2x[2] - x[1]*d3x[1] + x[2]*d3x[2]) + 32*x[1]*dx[3]*d3x[3] + x[1]^2*d4x[1] + 11*x[2]^2*d4x[1] - 10*x[1]*x[2]*d4x[2]) + 4*x[3]^2*(6*d2x[1]*d2x[3] + 8*dx[3]*d3x[1] - x[1]*d4x[3]) + x[2]^2*(-30*d2x[1]*d2x[3] - 32*dx[3]*d3x[1] - 8*dx[1]*d3x[3] + x[1]*d4x[3]) + x[1]*(96*dx[2]*dx[3]*d2x[2] + 12*dx[1]^2*d2x[3] + 36*dx[2]^2*d2x[3] + 6*(-8*dx[3]^2 + x[1]*d2x[1])*d2x[3] + 8*x[1]*dx[1]*d3x[3] + x[1]^2*d4x[3]))/15.

@inline D3_S333(x::MVector{3, Float64}, dx::MVector{3, Float64}, d2x::MVector{3, Float64}, d3x::MVector{3, Float64}, d4x::MVector{3, Float64})::Float64 = (6*(x[2]*d2x[1] - x[1]*d2x[2])*(dx[1]^2 + dx[2]^2 - 4*dx[3]^2 + x[1]*d2x[1] + x[2]*d2x[2] - 4*x[3]*d2x[3]) - 6*(x[1]*dx[1] + x[2]*dx[2] - 4*x[3]*dx[3])*(-(dx[2]*d2x[1]) + dx[1]*d2x[2] - x[2]*d3x[1] + x[1]*d3x[2]) + 2*(x[2]*dx[1] - x[1]*dx[2])*(3*dx[1]*d2x[1] + 3*dx[2]*d2x[2] - 12*dx[3]*d2x[3] + x[1]*d3x[1] + x[2]*d3x[2] - 4*x[3]*d3x[3]) - (x[1]^2 + x[2]^2 - 4*x[3]^2)*(-2*dx[2]*d3x[1] + 2*dx[1]*d3x[2] - x[2]*d4x[1] + x[1]*d4x[2]))/5.


end