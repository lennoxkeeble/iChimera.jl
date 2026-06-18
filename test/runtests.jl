using iChimera
using Test
using StaticArrays
using MAT

###### TEST HIGH-ORDER DERIVATIVES WRT BL TIME ######
dx_dt =  zeros(3);
d2x_dt =  zeros(3);
d3x_dt =  zeros(3);
d4x_dt =  zeros(3);
d5x_dt =  zeros(3);
d6x_dt =  zeros(3);
d7x_dt =  zeros(3);
d8x_dt =  zeros(3);
d9x_dt =  zeros(3);
dx_dλ = zeros(3);
d2x_dλ = zeros(3);
d3x_dλ = zeros(3);
d4x_dλ = zeros(3);
d5x_dλ = zeros(3);
d6x_dλ = zeros(3);
d7x_dλ = zeros(3);
d8x_dλ = zeros(3);
d9x_dλ = zeros(3);

a = 0.98; p = 8.0; e = 0.8; θmin = π /4; sign_Lz = 1;
x = [6.0, π / 3, π / 6];
E, L, Q, C = iChimera.ConstantsOfMotion.compute_ELC(a, p, e, θmin, sign_Lz) 

sign_dr = 1.0;
sign_dθ = 1.0;

dt_dλ, d2t_dλ, d3t_dλ, d4t_dλ, d5t_dλ, d6t_dλ, d7t_dλ, d8t_dλ, d9t_dλ, dλ_dt, d2λ_dt, d3λ_dt, d4λ_dt, d5λ_dt, d6λ_dt, d7λ_dt, d8λ_dt, d9λ_dt = iChimera.CoordinateDerivs.ComputeDerivs!(x, sign_dr, sign_dθ, dx_dt, d2x_dt, d3x_dt, d4x_dt, d5x_dt, d6x_dt, d7x_dt, d8x_dt, d9x_dt, dx_dλ, d2x_dλ, d3x_dλ, d4x_dλ, d5x_dλ, d6x_dλ, d7x_dλ, d8x_dλ, d9x_dλ, a, E, L, C);

dx_dt_MMA=[0.1492150850407141,0.037210239802416833,0.0678382150512224];
d2x_dt_MMA=[0.0054457060117013645,0.0001656341886566554,-0.005430468260244029];
d3x_dt_MMA=[-0.00033221256748209156,-0.0003664406381117368,0.0005908213533677006];
d4x_dt_MMA=[-5.320296612244223e-6,0.00008318009866873118,-0.00004410428427619115];
d5x_dt_MMA=[3.6464500725384553e-6,-0.000013085845529107754,-0.000012817504881439];
d6x_dt_MMA=[-1.0961644542936293e-7,3.530003021771099e-7,0.00001063457902129713];
d7x_dt_MMA=[-1.531612083481499e-7,1.0921908120100679e-6,-4.901596726845045e-6];
d8x_dt_MMA=[3.528463132436331e-8,-7.706478945590693e-7,1.5187797953974828e-6];
d9x_dt_MMA=[5.908769815132288e-9,3.6874883852653403e-7,8.341360851487352e-9];

@test isapprox(dx_dt, dx_dt_MMA, rtol = 1e-10)
@test isapprox(d2x_dt, d2x_dt_MMA, rtol = 1e-10)
@test isapprox(d3x_dt, d3x_dt_MMA, rtol = 1e-10)
@test isapprox(d4x_dt, d4x_dt_MMA, rtol = 1e-10)
@test isapprox(d5x_dt, d5x_dt_MMA, rtol = 1e-10)
@test isapprox(d6x_dt, d6x_dt_MMA, rtol = 1e-10)
@test isapprox(d7x_dt, d7x_dt_MMA, rtol = 1e-10)
@test isapprox(d8x_dt, d8x_dt_MMA, rtol = 1e-10)
@test isapprox(d9x_dt, d9x_dt_MMA, rtol = 1e-10)

##### TEST CONSISTENCY WITH HARMONIC COORDINATE TRANSFORMATIONS ######
xH = zeros(3);
dxH_dt = zeros(3);
d2xH_dt = zeros(3);
d3xH_dt = zeros(3);
d4xH_dt = zeros(3);
d5xH_dt = zeros(3);
d6xH_dt = zeros(3);
d7xH_dt = zeros(3);
d8xH_dt = zeros(3);
d9xH_dt = zeros(3);

iChimera.HarmonicCoordDerivs.compute_harmonic_derivs!(x, dx_dt, d2x_dt, d3x_dt, d4x_dt, d5x_dt, d6x_dt, d7x_dt, d8x_dt, d9x_dt, xH, dxH_dt, d2xH_dt, d3xH_dt, d4xH_dt, d5xH_dt, d6xH_dt, d7xH_dt, d8xH_dt, d9xH_dt, a);

xH ≈ iChimera.HarmonicCoords.xBLtoH(x, a)
dxH_dt ≈ iChimera.HarmonicCoords.vBLtoH(xH, dx_dt, a)
d2xH_dt ≈ iChimera.HarmonicCoords.aBLtoH(xH, dx_dt, d2x_dt, a)

d2x_dt ≈ iChimera.HarmonicCoords.aHtoBL(xH, dxH_dt, d2xH_dt, a)

iChimera.HarmonicCoords.aHtoBL(xH, zeros(3), d2xH_dt, a) ≈ iChimera.HarmonicCoords.jBLH(xH, a) * d2xH_dt

###### TESTING KERR METRIC AND KERR POTENTIALS in CARTESIAN COORDINATES ######
r, θ, ϕ = x;
r_dot, θ_dot, ϕ_dot = dx_dt;
xBL = x;
vBL = dx_dt;

xC = zero(xBL);
vC = zero(xBL);
iChimera.CartesianCoords.xBLtoH!(xC, xBL, a);
iChimera.CartesianCoords.vBLtoH!(vC, xC, vBL, a)

jBLC = iChimera.CartesianCoords.jBLH(xC, a)
jCBL = iChimera.CartesianCoords.jHBL(xC, a)
jBLC * jCBL ≈ [i == j ? 1 : 0 for i = 1:3, j = 1:3]   # identity_m

g_μν_cartesian = [-0.6688750858855246 0.023418944641674547 -0.04056280197902323 0;0.023418944641674547 1.2637721759916598 0.1330614203406636 0.16696288919982583;-0.04056280197902323 0.1330614203406636 1.1101260822867878 0.0963960690241971;0 0.16696288919982583 0.0963960690241971 1.1179780372443284];

@test isapprox(g_μν_cartesian, iChimera.CartesianCoords.g_μν_H(xC, a))

K_μν_cartesian = [0.3311249141144754 0.023418944641674547 -0.04056280197902323 0;0.023418944641674547 0.2637721759916598 0.1330614203406636 0.16696288919982583;-0.04056280197902323 0.1330614203406636 0.11012608228678766 0.0963960690241971;0 0.16696288919982583 0.0963960690241971 0.11797803724432843];

@test isapprox(K_μν_cartesian, iChimera.SelfAccelerationCartesian.K_μν(xC, a))
∂K_∂xk_cartesian = [-0.04029372032804981,-0.023263590278051048,-0.028325005467046065];
∂Ki_∂xk_cartesian = [-0.008704525969479886 0.00606272523960341 0; 0.003988395583691195 0.00870452596947989 0; -0.005906455007850681 0.010230280166217015 0]

∂Kij_∂x_cartesian = [-0.0002210978338527958 -0.024743268484635452 -0.032673280393962195;-0.024743268484635452 -0.039937185151102894 -0.04028527590280598;-0.032673280393962195 -0.04028527590280598 -0.049852151995948946]
∂Kij_∂y_cartesian = [-0.05926605993418717 0.019858043658625084 -0.04028527590280598;0.019858043658625084 0.036080664442629834 0.013844149379764591;-0.04028527590280598 0.013844149379764591 -0.028782153374543262]
∂Kij_∂z_cartesian = [-0.06880797430599606 -0.03555699542453498 0.009874314732905037;-0.03555699542453498 -0.02775029254613693 0.005700938269105817;0.009874314732905037 0.005700938269105817 0.045908888903379574]

rC = sqrt(xC[1]^2 + xC[2]^2 + xC[3]^2);
HessBLH = [iChimera.CartesianCoords.HessBLH(xC, rC, a, m) for m=1:3]
∂K_∂xk = @SVector [iChimera.SelfAccelerationCartesian.∂K_∂xk(xC, xBL, jBLC, HessBLH, a, j) for j=1:3];
∂Ki_∂xk = @SMatrix [iChimera.SelfAccelerationCartesian.∂Ki_∂xk(xC, rC, xBL, jBLC, HessBLH, a, j, k) for j=1:3, k=1:3];
∂Kij_∂xk = @SArray [iChimera.SelfAccelerationCartesian.∂Kij_∂xk(xC, rC, xBL, jBLC, HessBLH, a, j, k, l) for j=1:3, k=1:3, l=1:3]

@test isapprox(∂K_∂xk_cartesian, ∂K_∂xk)
@test isapprox(∂Ki_∂xk_cartesian, ∂Ki_∂xk)
@test isapprox(∂Kij_∂x_cartesian, ∂Kij_∂xk[1, :, :])
@test isapprox(∂Kij_∂y_cartesian, ∂Kij_∂xk[2, :, :])
@test isapprox(∂Kij_∂z_cartesian, ∂Kij_∂xk[3, :, :])

###### TESTING KERR METRIC AND KERR POTENTIALS in HARMONIC COORDINATES ######
xH = iChimera.HarmonicCoords.xBLtoH(xBL, a);
vH = iChimera.HarmonicCoords.vBLtoH(xH, vBL, a);

@test iChimera.HarmonicCoords.xHtoBL(xH, a) ≈ xBL
@test iChimera.HarmonicCoords.vHtoBL(xH, vH, a) ≈ vBL

jBLH = iChimera.HarmonicCoords.jBLH(xH, a)
jHBL = iChimera.HarmonicCoords.jHBL(xH, a)
@test jBLH * jHBL ≈ [i == j ? 1 : 0 for i = 1:3, j = 1:3]   # identity_m

g_μν_harmonic = [-0.6688750858855246 0.02742183732650591 -0.047857416118385475 -0.00018929217690248146; 0.02742183732650591 1.4598127412586643 0.02398293333694035 0.019295596813892858; -0.047857416118385475 0.02398293333694035 1.4516925363102926 0.01689304425420316; -0.00018929217690248146 0.019295596813892858 0.01689304425420316 1.450622640701141];

@test isapprox(g_μν_harmonic, iChimera.HarmonicCoords.g_μν_H(xH, a))

K_μν_harmonic = [0.3311249141144754 0.02742183732650591 -0.047857416118385475 -0.00018929217690248146; 0.02742183732650591 0.45981274125866456 0.02398293333694035 0.019295596813892858; -0.047857416118385475 0.02398293333694035 0.4516925363102924 0.01689304425420316; -0.00018929217690248146 0.019295596813892858 0.01689304425420316 0.45062264070114116];

@test isapprox(K_μν_harmonic, iChimera.SelfAccelerationHarmonic.K_μν(xH, a))

∂K_∂xk_harmonic = [-0.04049907038065751, -0.023520435081995278, -0.029308896742549204];
∂Ki_∂xk_harmonic = [-0.01131445893792844 0.007344044731051753 0.00014615832561827783; 0.005971207711286233 0.011451749399193707 0.00008488361281101638; -0.007717304730092342 0.013798959435045813 0.00007664485119975283];

∂Kij_∂x_harmonic = [-0.06854260928760211 -0.007703297461706768 -0.0064355247140344685; -0.007703297461706768 -0.07932674739904806 -0.010474326770334053; -0.0064355247140344685 -0.010474326770334053 -0.07913319895492536];
∂Kij_∂y_harmonic = [-0.052377805309502756 -0.002345698857464079 -0.008164776190600187; -0.002345698857464079 -0.03349951414190476 -0.0010262208473323375; -0.008164776190600187 -0.0010262208473323375 -0.04595777758244239];
∂Kij_∂z_harmonic = [-0.06031240159552314 -0.011255322845189597 -0.0007406869120335224; -0.011255322845189597 -0.0579530065454595 -0.0012964347788881753; -0.0007406869120335224 -0.0012964347788881753 -0.045460456754370826];


rH = sqrt(xH[1]^2 + xH[2]^2 + xH[3]^2);
HessBLH = [iChimera.HarmonicCoords.HessBLH(xH, rH, a, m) for m=1:3]
∂K_∂xk = @SVector [iChimera.SelfAccelerationHarmonic.∂K_∂xk(xH, xBL, jBLH, HessBLH, a, j) for j=1:3];
∂Ki_∂xk = @SMatrix [iChimera.SelfAccelerationHarmonic.∂Ki_∂xk(xH, rH, xBL, jBLH, HessBLH, a, j, k) for j=1:3, k=1:3];
∂Kij_∂xk = @SArray [iChimera.SelfAccelerationHarmonic.∂Kij_∂xk(xH, rH, xBL, jBLH, HessBLH, a, j, k, l) for j=1:3, k=1:3, l=1:3]

@test isapprox(∂K_∂xk_harmonic, ∂K_∂xk)
@test isapprox(∂Ki_∂xk_harmonic, ∂Ki_∂xk)
@test isapprox(∂Kij_∂x_harmonic, ∂Kij_∂xk[1, :, :])
@test isapprox(∂Kij_∂y_harmonic, ∂Kij_∂xk[2, :, :])
@test isapprox(∂Kij_∂z_harmonic, ∂Kij_∂xk[3, :, :])

###### TESTING GAMMA EQUAL IN HARMONIC AND BL COORDINATES ######
@test iChimera.SelfAccelerationHarmonic.Γ(vH, xH, a) ≈ iChimera.BLTimeGeodesics.Γ(r, θ, ϕ, vBL, a)
@test iChimera.SelfAccelerationCartesian.Γ(vC, xC, a) ≈ iChimera.BLTimeGeodesics.Γ(r, θ, ϕ, vBL, a)


###### TESTING COORDINATE DERIVATIVE COMPUTATION IN HARMONIC COORDINATES ######
xHDerivs_MMA = [1.633536310632804, 3.6207825800724716, 2.6590173841405575, -25.477583917552835, -198.22937662413975, -604.3710144278622, -2103.685204284686, 50196.725027620996, -302035.7873468707, 1.3420726353827521*1e7];
yHDerivs_MMA = [0.7204575781016996, 2.6037850458877645, 9.433240825620512, 14.424221724291963, -25.315039292192818, -1100.1181048031633, -2817.6511045173906, -78425.90005093698, 786760.7841776566, -1.5944788233137697*1e7];
zHDerivs_MMA = [1., -0.6411496334238431, -11.243389142737914, -43.77408426449627, -106.78248444718885, 11.278780437300327, 2326.929754644706, 19606.415106847817, 113533.89599864538,  464196.5937074214];

begin
    rD(n) = 3 + n/2
    thetaD(n) = π/3 + n/3
    phiD(n) = π/8 + n/4

    xBL = [rD(0), thetaD(0), phiD(0)];
    dxBL = [rD(1), thetaD(1), phiD(1)];
    d2xBL = [rD(2), thetaD(2), phiD(2)];
    d3xBL = [rD(3), thetaD(3), phiD(3)];
    d4xBL = [rD(4), thetaD(4), phiD(4)];
    d5xBL = [rD(5), thetaD(5), phiD(5)];
    d6xBL = [rD(6), thetaD(6), phiD(6)];
    d7xBL = [rD(7), thetaD(7), phiD(7)];
    d8xBL = [rD(8), thetaD(8), phiD(8)];
    d9xBL = [rD(9), thetaD(9), phiD(9)];

    xH = zeros(3);
    dxH = zeros(3);
    d2xH = zeros(3);
    d3xH = zeros(3);
    d4xH = zeros(3);
    d5xH = zeros(3);
    d6xH = zeros(3);
    d7xH = zeros(3);
    d8xH = zeros(3);
    d9xH = zeros(3);
    a = 0.5;

    xH_MMA = [xHDerivs_MMA[1], yHDerivs_MMA[1], zHDerivs_MMA[1]]
    dxH_MMA = [xHDerivs_MMA[2], yHDerivs_MMA[2], zHDerivs_MMA[2]]
    d2xH_MMA = [xHDerivs_MMA[3], yHDerivs_MMA[3], zHDerivs_MMA[3]]
    d3xH_MMA = [xHDerivs_MMA[4], yHDerivs_MMA[4], zHDerivs_MMA[4]]
    d4xH_MMA = [xHDerivs_MMA[5], yHDerivs_MMA[5], zHDerivs_MMA[5]]
    d5xH_MMA = [xHDerivs_MMA[6], yHDerivs_MMA[6], zHDerivs_MMA[6]]
    d6xH_MMA = [xHDerivs_MMA[7], yHDerivs_MMA[7], zHDerivs_MMA[7]]
    d7xH_MMA = [xHDerivs_MMA[8], yHDerivs_MMA[8], zHDerivs_MMA[8]]
    d8xH_MMA = [xHDerivs_MMA[9], yHDerivs_MMA[9], zHDerivs_MMA[9]]
    d9xH_MMA = [xHDerivs_MMA[10], yHDerivs_MMA[10], zHDerivs_MMA[10]]

    iChimera.HarmonicCoordDerivs.compute_harmonic_derivs!(xBL, dxBL, d2xBL, d3xBL, d4xBL, d5xBL, d6xBL, d7xBL, d8xBL, d9xBL, xH, dxH, d2xH, d3xH, d4xH, d5xH, d6xH, d7xH, d8xH, d9xH, a)

    @test isapprox(dxH_MMA, dxH) && isapprox(d2xH_MMA, d2xH) && isapprox(d3xH_MMA, d3xH) && isapprox(d4xH_MMA, d4xH) && isapprox(d5xH_MMA, d5xH) && isapprox(d6xH_MMA, d6xH) && isapprox(d7xH_MMA, d7xH) && isapprox(d8xH_MMA, d8xH) && isapprox(d9xH_MMA, d9xH)
end

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

v2power = 0;
Vrr, ∂Vrr_∂t = iChimera.RRPotentials.compute_RR_potentials!(Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a, x, dx, d2x, Mij5, Mij6, Mij7, Mij8, dxmMij5, dxmMij6, dxmMij7, Mijk7, Mijk8, dxmMijk7, Sij5, Sij6, dxmSij5, v2power; rr_derivative_model=:legacy_worldline);

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

###### TESTING LOW ORDER TIME DERIVATIVES OF MULTIPOLE MOMENTS DONE BY HAND ######
δ(i, j) = i == j ? 1 : 0
Mij2(x, dx, d2x, i::Int64, j::Int64)::Float64 = -(2.0 / 3.0) * δ(i, j) * (x[1] * d2x[1] + dx[1]^2 + x[2] * d2x[2] + dx[2]^2 + x[3] * d2x[3] + dx[3]^2) + x[j] * d2x[i] + 2 * dx[i] * dx[j] + x[i] * d2x[j]

Mij2_test = [Mij2(x, 2 * dx, 3* d2x, i, j) for i=1:3, j=1:3]

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

aSF_H = zeros(4)
aSF_BL = zeros(4)
iChimera.SelfAccelerationCartesian.aRRα(aSF_H, aSF_BL, xC, v, vC, xBL, rC, a, Vrr, ∂Vrr_∂t, Virr, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a)
@test -iChimera.SelfAccelerationCartesian.Γ(vC, xC, a)^2 * iChimera.SelfAccelerationCartesian.Pαβ(vC, xC, a) * (iChimera.SelfAccelerationCartesian.A1_β(v, vC, ∂Vrr_∂t, ∂Vrr_∂a, ∂Virr_∂t, ∂Virr_∂a) + iChimera.SelfAccelerationCartesian.A2_β(xC, vC, xBL, rC, a, Vrr, Virr)) ≈ aSF_H


###### test harmonic coordinates and their derivatives in the far-field limit where they should (NOT) reduce to flat-space Cartesian coordinates ######
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

include("test_harmonic_rr_acceleration.jl")
include("test_rr_derivative_model.jl")


include("/Users/lennoxkeeble/KerrMetricDedoner.jl")

xH = [5.0, 10.0, 15.0]; a = 0.5;
g_μν = iChimera.HarmonicCoords.g_μν_H(xH, a)
gμν = iChimera.HarmonicCoords.gμν_H(xH, a)

m = 1.0
# they use the + - - -  signature
g_μν_2 = -KerrMetricDeDonder.get_metric_components(xH..., m, a)
gμν_2 = -KerrMetricDeDonder.get_inverse_metric_components(xH..., m, a)

g_μν ≈ g_μν_2
gμν ≈ gμν_2

gμν - gμν_2

g_μν * gμν

gμν_2 * g_μν_2
