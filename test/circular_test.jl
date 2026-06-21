using iChimera
using LinearAlgebra
using StaticArrays
using Test

const HRR = iChimera.HarmonicRRAcceleration
const KM = iChimera.Kerr.KerrMetric

function unconstrained_coordinate_time_fluxes(gBL, a4_BL, u_t)
    a_t = gBL[1, 1] * a4_BL[1] + gBL[1, 4] * a4_BL[4]
    a_phi = gBL[4, 1] * a4_BL[1] + gBL[4, 4] * a4_BL[4]
    Edot = -a_t / u_t
    Ldot = a_phi / u_t
    return Edot, Ldot
end

function circular_energy_angular_momentum_balance_lhs(;
    spin::Float64=0.7,
    p::Float64=10.0,
    phi::Float64=0.3,
    q::Float64=1.0e-5,
)
    theta = pi / 2
    omega_phi = 1.0 / (spin + p * sqrt(p))

    xBL = @SVector [p, theta, phi]
    vBL = @SVector [0.0, 0.0, omega_phi]
    xH = iChimera.HarmonicCoords.xBLtoH(xBL, spin)
    vH = iChimera.HarmonicCoords.vBLtoH(xH, vBL, spin)

    Acoord_H = HRR.rr_harmonic_accel(xH, vH, q; include25=true, include35=true)
    gH = iChimera.HarmonicCoords.g_μν_H(xH, spin)
    a4_H = HRR.lift_coord_accel_to_four(Acoord_H, vH, gH)
    a4_BL = HRR.four_accel_H_to_BL(a4_H, xH, spin)

    gBL = @SMatrix [KM.g_μν(p, theta, phi, spin, mu, nu) for mu in 1:4, nu in 1:4]
    u_t = HRR.gamma_from_metric_velocity(gBL, vBL)
    u4_BL = @SVector [u_t, 0.0, 0.0, u_t * omega_phi]
    u4_BL_cov = gBL * u4_BL

    # Do not call Evolve_BL here: its circular-equatorial branch imposes
    # Ldot = Edot / omega_phi. This is the raw numerical check using only
    # Edot = -a_t / u^t and Ldot = a_phi / u^t.
    Edot, Ldot = unconstrained_coordinate_time_fluxes(gBL, a4_BL, u_t)

    lhs = Edot - omega_phi * Ldot
    rhs = -dot(u4_BL_cov, a4_BL) / u_t^2

    return (
        lhs = lhs,
        rhs = rhs,
        u_dot_a = dot(u4_BL_cov, a4_BL),
        Edot = Edot,
        Ldot = Ldot,
        omega_phi = omega_phi,
    )
end

@testset "Circular-orbit energy-angular-momentum balance" begin
    result = circular_energy_angular_momentum_balance_lhs()

    println("Circular-orbit Edot - Ωphi Ldot = $(result.lhs)")
    println("Circular-orbit -u_alpha a^alpha / (u^t)^2 = $(result.rhs)")
    println("Circular-orbit u_alpha a^alpha = $(result.u_dot_a)")

    @test isapprox(result.lhs, result.rhs; atol=1.0e-20, rtol=1.0e-12)
    @test isapprox(result.lhs, 0.0; atol=1.0e-20)
end
