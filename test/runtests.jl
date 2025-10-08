using iChimera
using Test
using StaticArrays
using MAT

###### TEST RADIATION REACTION POTENTIALS ######
DX1(n) = sin(n+1) * cos(n+1)
DX2(n) = sin(5n+2) * cos(5n+2)
DX3(n) = sin(10n+3) * cos(10n+3)

x = @MArray [DX1(0), DX2(0), DX3(0)];
dx = @MArray [DX1(1), DX2(1), DX3(1)];
d2x = @MArray [DX1(2), DX2(2), DX3(2)];
d3x = @MArray [DX1(3), DX2(3), DX3(3)];
d4x = @MArray [DX1(4), DX2(4), DX3(4)];
d5x = @MArray [DX1(5), DX2(5), DX3(5)];
d6x = @MArray [DX1(6), DX2(6), DX3(6)];
d7x = @MArray [DX1(7), DX2(7), DX3(7)];
d8x = @MArray [DX1(8), DX2(8), DX3(8)];
d9x = @MArray [DX1(9), DX2(9), DX3(9)];

Virr = @MArray zeros(3)
∂Vrr_∂a = @MArray zeros(3)
∂Virr_∂t = @MArray zeros(3)
∂Virr_∂a = @MMatrix zeros(3, 3)
q = 0.0042839;

# note that multipole moments and their derivatives already checked in its own package
Mij5 = @MArray zeros(3, 3)
Mij6 = @MArray zeros(3, 3)
Mij7 = @MArray zeros(3, 3)
Mij8 = @MArray zeros(3, 3)
dxmMij5 = @MArray zeros(3, 3, 3)
dxmMij6 = @MArray zeros(3, 3, 3)
dxmMij7 = @MArray zeros(3, 3, 3)

Mijk7 = @MArray zeros(3, 3, 3)
Mijk8 = @MArray zeros(3, 3, 3)
dxmMijk7 = @MArray zeros(3, 3, 3, 3)

Sij5 = @MArray zeros(3, 3)
Sij6 = @MArray zeros(3, 3)
dxmSij5 = @MArray zeros(3, 3, 3)

OnePN, TwoPN, TwoPointFivePN = 1.0, 1.0, 1.0;
iChimera.MultipoleDerivs.compute_SF_moments!(q, Mij5, Mij6, Mij7, Mij8, dxmMij5, dxmMij6, dxmMij7, Mijk7, Mijk8, dxmMijk7, Sij5, Sij6, dxmSij5, x, dx, d2x, d3x, d4x, d5x, d6x, d7x, d8x, d9x, OnePN, TwoPN, TwoPointFivePN)

Virr = @MVector zeros(3)
∂Vrr_∂a = @MVector zeros(3)
∂Virr_∂t = @MVector zeros(3)
∂Virr_∂a = @MArray zeros(3, 3)

Vrr, ∂Vrr_∂t = iChimera.RRPotentials.compute_RR_potentials!(Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a, x, dx, Mij5, Mij6, Mij7, Mij8, dxmMij5, dxmMij6, dxmMij7, Mijk7, Mijk8, dxmMijk7, Sij5, Sij6, dxmSij5);

# load results from explicit MMA computation
function get_mata_arr(fname::String)
    file = matopen("./test/MMA/"*fname*".mat")
    arr_julia = read(file, "Expression1")  # variable name must match
    close(file)
    return arr_julia
end

Vrr_MMA = get_mata_arr("Vrr");
dtVrr_MMA = get_mata_arr("dtVrr");
dxVrr_MMA = get_mata_arr("dxVrr");
Virr_MMA = get_mata_arr("Virr");
dtVirr_MMA = get_mata_arr("dtVirr");
dxVirr_MMA = get_mata_arr("dxVirr");

@test Vrr ≈ Vrr_MMA
@test ∂Vrr_∂t ≈ dtVrr_MMA
@test ∂Vrr_∂a ≈ dxVrr_MMA
@test Virr ≈ Virr_MMA
@test ∂Virr_∂t ≈ dtVirr_MMA
@test ∂Virr_∂a ≈ dxVirr_MMA

###### TEST SELF-ACCELERATION COMPUTATION ######
# TEST A1 piece
vH = dx;
v = sqrt(vH[1]^2 + vH[2]^2 + vH[3]^2);
∂Vrr_∂t = 0.1;
∂Vrr_∂a = @MArray [0.2, -0.3, 0.4];
∂Virr_∂t = @MArray [0.05, -0.06, 0.07];
∂Virr_∂a = @MMatrix [0.01 0.02 0.03; -0.04 0.05 -0.06; 0.07 -0.08 0.09];

A1 = iChimera.SelfAccelerationHarmonic.A1_β(v, vH, ∂Vrr_∂t, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a)

A1_MMA = @MArray [-0.05556309318452976,
-0.5860179813375077,
0.6068225738467036,
-0.7718514324089021]
@test isapprox(A1, A1_MMA, rtol=1e-14)

# content of function A2_β
∂K_∂xk = @SVector [0.5, 0.6, 0.7];
∂Ki_∂xk = @SMatrix [2.0 * (i + j) for i=1:3, j=1:3];
∂Kij_∂xk = @SArray [3.0 * (j + k + l) for j=1:3, k=1:3, l=1:3];
Q = 1.1;
Qi = @SVector [1.2, 1.3, 1.4];
Qij = @SMatrix [4.0 * (i + j) for i=1:3, j=1:3];

BRR = iChimera.SelfAccelerationHarmonic.B_RR(Qi, ∂K_∂xk)
BiRR = iChimera.SelfAccelerationHarmonic.Bi_RR(Qij, ∂K_∂xk)

CRR = iChimera.SelfAccelerationHarmonic.C_RR(vH, ∂K_∂xk, ∂Ki_∂xk, Q, Qi)
CiRR = iChimera.SelfAccelerationHarmonic.Ci_RR(vH, ∂K_∂xk, ∂Ki_∂xk, Qi, Qij)

DRR = iChimera.SelfAccelerationHarmonic.D_RR(vH, ∂Ki_∂xk, ∂Kij_∂xk, Q, Qi)
DiRR = iChimera.SelfAccelerationHarmonic.Di_RR(vH, ∂Ki_∂xk, ∂Kij_∂xk, Qi, Qij)

BRR_MMA = 2.36;
BiRR_MMA = @SVector [-45.8, -60.4, -75.0];
CRR_MMA = -0.07497540810986807;
CiRR_MMA = @SVector [1.799409794636832, 1.9493606108565682, 2.099311427076304];
DRR_MMA = -27.127296024255514;
DiRR_MMA = @SVector [525.8385912324655, 690.441408668543, 855.0442261046204];


@test isapprox(BRR, BRR_MMA, rtol=1e-14)
@test isapprox(BiRR, BiRR_MMA, rtol=1e-14)
@test isapprox(CRR, CRR_MMA, rtol=1e-14)
@test isapprox(CiRR, CiRR_MMA, rtol=1e-14)
@test isapprox(DRR, DRR_MMA, rtol=1e-14)
@test isapprox(DiRR, DiRR_MMA, rtol=1e-14)

VRR = π/3;
ViRR = @SVector [π/6, -π/5, π/4];

A2_t = (BRR + CRR + DRR) * VRR + iChimera.SelfAccelerationHarmonic.dot((BiRR + CiRR + DiRR), ViRR)   # Eq. 62
A2_i = -2.0 * (BRR + CRR + DRR) * ViRR - (BiRR + CiRR + DiRR) * VRR / 2.0  # Eq. 63

A2_MMA = @SVector [443.47760787402143, -226.2750215643884, -362.1273119562816, -370.50724987851686];

@test isapprox(vcat(A2_t, A2_i), A2_MMA, rtol=1e-14)

###### test harmonic coordinates and their derivatives in the far-field limit where they should reduce to flat-space Cartesian coordinates ######
xBL = @MArray zeros(3); vBL = @MArray zeros(3); aBL = @MArray zeros(3);

dxBL_dt = @MArray zeros(3); d2xBL_dt = @MArray zeros(3); d3xBL_dt = @MArray zeros(3); d4xBL_dt = @MArray zeros(3);
d5xBL_dt = @MArray zeros(3); d6xBL_dt = @MArray zeros(3); d7xBL_dt = @MArray zeros(3); d8xBL_dt = @MArray zeros(3); d9xBL_dt = @MArray zeros(3);

dx_dλ = @MArray zeros(3); d2x_dλ = @MArray zeros(3); d3x_dλ = @MArray zeros(3); d4x_dλ = @MArray zeros(3);
d5x_dλ = @MArray zeros(3); d6x_dλ = @MArray zeros(3); d7x_dλ = @MArray zeros(3); d8x_dλ = @MArray zeros(3); d9x_dλ = @MArray zeros(3);

xH = @MArray zeros(3); dxH_dt = @MArray zeros(3); d2xH_dt = @MArray zeros(3); d3xH_dt = @MArray zeros(3); d4xH_dt = @MArray zeros(3);
d5xH_dt = @MArray zeros(3); d6xH_dt = @MArray zeros(3); d7xH_dt = @MArray zeros(3); d8xH_dt = @MArray zeros(3); d9xH_dt = @MArray zeros(3);

xC = @MArray zeros(3); dxC_dt = @MArray zeros(3); d2xC_dt = @MArray zeros(3); d3xC_dt = @MArray zeros(3); d4xC_dt = @MArray zeros(3);
d5xC_dt = @MArray zeros(3); d6xC_dt = @MArray zeros(3); d7xC_dt = @MArray zeros(3); d8xC_dt = @MArray zeros(3); d9xC_dt = @MArray zeros(3);

vH = @MArray zeros(3);
aH = @MArray zeros(3);

a = 0.5; e = 0.4; θmin = π/6; sign_Lz = 1;
p_vals = [1.0e1, 1.0e2, 500.0]
print_errors = true

for p in p_vals

    # calculate integrals of motion from orbital parameters
    EEi, LLi, CCi = iChimera.ConstantsOfMotion.SchmidtELQ(a, p, e, θmin, sign_Lz)
    QQi = CCi + (LLi - a * EEi)^2

    rplus = iChimera.Kerr.KerrMetric.rplus(a); rminus = iChimera.Kerr.KerrMetric.rminus(a);

    # periastron and apastron
    rp = p / (1 + e)
    ra = p / (1 - e)

    # compute roots of radial function R(r)
    zm = cos(θmin)^2
    zp = CCi / (a^2 * (1.0-EEi^2) * zm)    # Eq. E23
    A = 1.0 / (1.0 - EEi^2) - (ra + rp) / 2.0    # Eq. E20
    B = a^2 * CCi / ((1.0 - EEi^2) * ra * rp)    # Eq. E21
    r3 = A + sqrt(A^2 - B); r4 = A - sqrt(A^2 - B);    # Eq. E19
    p3 = r3 * (1.0 - e); p4 = r4 * (1.0 + e);  # Above Eq. 96

    rr = ra; θθ = θmin; ϕϕ = π/3;
    psi = 1.0 * π; chi = 0.0

    # compute time derivatives
    psi_dot = iChimera.BLTimeGeodesics.psi_dot(psi, chi, ϕϕ, a, EEi, LLi, p, e, θmin, p3, p4, zp, zm)
    chi_dot = iChimera.BLTimeGeodesics.chi_dot(psi, chi, ϕϕ, a, EEi, LLi, p, e, θmin, p3, p4, zp, zm)
    ϕ_dot = iChimera.BLTimeGeodesics.phi_dot(psi, chi, ϕϕ, a, EEi, LLi, p, e, θmin, p3, p4, zp, zm)

    # compute BL coordinates t, r, θ and their time derivatives
    r = iChimera.BLTimeGeodesics.r(psi, p, e)
    θ = acos((π/2<chi<1.5π) ? -sqrt(iChimera.BLTimeGeodesics.z(chi, θmin)) : sqrt(iChimera.BLTimeGeodesics.z(chi, θmin)))
    r ≈ rr
    θ ≈ θθ

    r_dot = iChimera.BLTimeGeodesics.dr_dt(psi_dot, psi, p, e);
    θ_dot = iChimera.BLTimeGeodesics.dθ_dt(chi_dot, chi, θ, θmin);
    v = [r_dot, θ_dot, ϕ_dot];
    dt_dτ = iChimera.BLTimeGeodesics.Γ(r, θ, ϕϕ, v, a)

    # substitute solution back into geodesic equation to find second derivatives of BL coordinates (wrt t)
    r_ddot = iChimera.BLTimeGeodesics.dr2_dt2(r, θ, ϕϕ, r_dot, θ_dot, ϕ_dot, a)
    θ_ddot = iChimera.BLTimeGeodesics.dθ2_dt2(r, θ, ϕϕ, r_dot, θ_dot, ϕ_dot, a)
    ϕ_ddot = iChimera.BLTimeGeodesics.dϕ2_dt2(r, θ, ϕϕ, r_dot, θ_dot, ϕ_dot, a)


    # COMPUTE BL COORDINATE DERIVATIVES
    xBL[1] = rr; xBL[2] = θθ; xBL[3] = ϕϕ;
    vBL[1] = r_dot; vBL[2] = θ_dot; vBL[3] = ϕ_dot;
    aBL[1] = r_ddot; aBL[2] = θ_ddot; aBL[3] = ϕ_ddot;
    iChimera.CoordinateDerivs.ComputeDerivs!(xBL, sign(vBL[1]), sign(vBL[2]), dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, d9xBL_dt, dx_dλ, d2x_dλ, d3x_dλ, d4x_dλ, d5x_dλ, d6x_dλ, d7x_dλ, d8x_dλ, d9x_dλ, a, EEi, LLi, CCi);

    # COMPUTE HARMONIC COORDINATE DERIVATIVES
    iChimera.HarmonicCoordDerivs.compute_harmonic_derivs!(xBL, dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, d9xBL_dt, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt, d5xH_dt, d6xH_dt, d7xH_dt, d8xH_dt, d9xH_dt, a);
    iChimera.CartesianCoords.compute_cartesian_derivs!(xBL, dxBL_dt, d2xBL_dt, d3xBL_dt, d4xBL_dt, d5xBL_dt, d6xBL_dt, d7xBL_dt, d8xBL_dt, d9xBL_dt, xC, dxC_dt, d2xC_dt, d3xC_dt, d4xC_dt, d5xC_dt, d6xC_dt, d7xC_dt, d8xC_dt, d9xC_dt);

    # println("r = $r, Δx = $(100 * (1 - xH[1] / xC[1]))%, Δy = $(100 * (1 - xH[2] / xC[2]))%, Δz = $(100 * (1 - xH[3] / xC[3]))%")
    # println("Δdx = $(100 * (1 - dxH_dt[1] / dxC_dt[1]))%, Δdy = $(100 * (1 - dxH_dt[2] / dxC_dt[2]))%, Δdz = $(100 * (1 - dxH_dt[3] / dxC_dt[3]))%, dxzH = $(dxH_dt[3]), dxzC = $(dxC_dt[3])")
    # println("Δd2x = $(100 * (1 - d2xH_dt[1] / d2xC_dt[1]))%, Δd2y = $(100 * (1 - d2xH_dt[2] / d2xC_dt[2]))%, Δd2z = $(100 * (1 - d2xH_dt[3] / d2xC_dt[3]))%, d2zH = $(d2xH_dt[3]), d2zC = $(d2xC_dt[3])")
    # println("Δd3x = $(100 * (1 - d3xH_dt[1] / d3xC_dt[1]))%, Δd3y = $(100 * (1 - d3xH_dt[2] / d3xC_dt[2]))%, Δd3z = $(100 * (1 - d3xH_dt[3] / d3xC_dt[3]))%, d3zH = $(d3xH_dt[3]), d3zC = $(d3xC_dt[3])")
    # println("Δd4x = $(100 * (1 - d4xH_dt[1] / d4xC_dt[1]))%, Δd4y = $(100 * (1 - d4xH_dt[2] / d4xC_dt[2]))%, Δd4z = $(100 * (1 - d4xH_dt[3] / d4xC_dt[3]))%, d4zH = $(d4xH_dt[3]), d4zC = $(d4xC_dt[3])")
    # println("Δd5x = $(100 * (1 - d5xH_dt[1] / d5xC_dt[1]))%, Δd5y = $(100 * (1 - d5xH_dt[2] / d5xC_dt[2]))%, Δd5z = $(100 * (1 - d5xH_dt[3] / d5xC_dt[3]))%, d5zH = $(d5xH_dt[3]), d5zC = $(d5xC_dt[3])")
    # println("Δd6x = $(100 * (1 - d6xH_dt[1] / d6xC_dt[1]))%, Δd6y = $(100 * (1 - d6xH_dt[2] / d6xC_dt[2]))%, Δd6z = $(100 * (1 - d6xH_dt[3] / d6xC_dt[3]))%, d6zH = $(d6xH_dt[3]), d6zC = $(d6xC_dt[3])")
    # println("Δd7x = $(100 * (1 - d7xH_dt[1] / d7xC_dt[1]))%, Δd7y = $(100 * (1 - d7xH_dt[2] / d7xC_dt[2]))%, Δd7z = $(100 * (1 - d7xH_dt[3] / d7xC_dt[3]))%, d7zH = $(d7xH_dt[3]), d7zC = $(d7xC_dt[3])")
    # println("Δd8x = $(100 * (1 - d8xH_dt[1] / d8xC_dt[1]))%, Δd8y = $(100 * (1 - d8xH_dt[2] / d8xC_dt[2]))%, Δd8z = $(100 * (1 - d8xH_dt[3] / d8xC_dt[3]))%, d8zH = $(d8xH_dt[3]), d8zC = $(d8xC_dt[3])\n")

    if print_errors
        println("r = $r, Δx = $(100 * (1 - xH[1] / xC[1]))%, Δy = $(100 * (1 - xH[2] / xC[2]))%, Δz = $(100 * (1 - xH[3] / xC[3]))%")
        println("Δdx = $(100 * (1 - dxH_dt[1] / dxC_dt[1]))%, Δdy = $(100 * (1 - dxH_dt[2] / dxC_dt[2]))%, Δdz = $(100 * (1 - dxH_dt[3] / dxC_dt[3]))%")
        println("Δd2x = $(100 * (1 - d2xH_dt[1] / d2xC_dt[1]))%, Δd2y = $(100 * (1 - d2xH_dt[2] / d2xC_dt[2]))%, Δd2z = $(100 * (1 - d2xH_dt[3] / d2xC_dt[3]))%")
        println("Δd3x = $(100 * (1 - d3xH_dt[1] / d3xC_dt[1]))%, Δd3y = $(100 * (1 - d3xH_dt[2] / d3xC_dt[2]))%, Δd3z = $(100 * (1 - d3xH_dt[3] / d3xC_dt[3]))%")
        println("Δd4x = $(100 * (1 - d4xH_dt[1] / d4xC_dt[1]))%, Δd4y = $(100 * (1 - d4xH_dt[2] / d4xC_dt[2]))%, Δd4z = $(100 * (1 - d4xH_dt[3] / d4xC_dt[3]))%")
        println("Δd5x = $(100 * (1 - d5xH_dt[1] / d5xC_dt[1]))%, Δd5y = $(100 * (1 - d5xH_dt[2] / d5xC_dt[2]))%, Δd5z = $(100 * (1 - d5xH_dt[3] / d5xC_dt[3]))%")
        println("Δd6x = $(100 * (1 - d6xH_dt[1] / d6xC_dt[1]))%, Δd6y = $(100 * (1 - d6xH_dt[2] / d6xC_dt[2]))%, Δd6z = $(100 * (1 - d6xH_dt[3] / d6xC_dt[3]))%")
        println("Δd7x = $(100 * (1 - d7xH_dt[1] / d7xC_dt[1]))%, Δd7y = $(100 * (1 - d7xH_dt[2] / d7xC_dt[2]))%, Δd7z = $(100 * (1 - d7xH_dt[3] / d7xC_dt[3]))%")
        println("Δd8x = $(100 * (1 - d8xH_dt[1] / d8xC_dt[1]))%, Δd8y = $(100 * (1 - d8xH_dt[2] / d8xC_dt[2]))%, Δd8z = $(100 * (1 - d8xH_dt[3] / d8xC_dt[3]))%\n")
    end
end

vH .= dxH_dt
aH .= d2xH_dt
rH = iChimera.SelfAccelerationHarmonic.norm_3d(xH);
v = iChimera.SelfAccelerationHarmonic.norm_3d(vH);
aSF_BL = @MVector zeros(4)
aSF_H = @MVector zeros(4)
iChimera.SelfAccelerationHarmonic.aRRα(aSF_H, aSF_BL, xH, v, vH, xBL, rH, a, Vrr, ∂Vrr_∂t, Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a)
iChimera.SelfAccelerationCartesian.aRRα(aSF_H, aSF_BL, xH, v, vH, xBL, rH, a, Vrr, ∂Vrr_∂t, Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a)

@testset "iChimera.jl" begin
    # Write your tests here.
end
