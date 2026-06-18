using iChimera
using Test
using StaticArrays
using LinearAlgebra

const HRR = iChimera.HarmonicRRAcceleration

@testset "Harmonic RR acceleration formulas" begin
    X = @SVector [8.1, -3.2, 4.7]
    V = @SVector [0.07, -0.03, 0.02]
    q = 1.0e-5

    A25 = HRR.rr_harmonic_accel(X, V, q; include25=true, include35=false)

    m = 1.0 * (1 + q)
    ν = q / (1 + q)^2
    R = norm(X)
    N = X / R
    Rdot = dot(N, V)
    V2 = dot(V, V)
    u = m / R

    pref25 = m^2 * ν / R^3
    expected25 = pref25 * (
        ((136 / 15) * u * Rdot + (24 / 5) * Rdot * V2) * N -
        ((24 / 5) * u + (8 / 5) * V2) * V
    )

    @test A25 ≈ expected25

    # This is the same 2.5PN harmonic branch written in the Iyer-Will
    # coefficient form. Different references swap gauge labels, so the test
    # checks the coefficient convention directly instead of relying on names.
    coeff_form = (8 / 5) * m^2 * ν / R^3 * (Rdot * (3 * V2 + (17 / 3) * u) * N - (V2 + 3 * u) * V)
    @test A25 ≈ coeff_form

    Xcirc = @SVector [3.0, 4.0, 1.0]
    Vcirc = 1.0e-2 .* @SVector [-4.0, 3.0, 0.0]
    @test isapprox(dot(Xcirc, Vcirc), 0.0; atol=1.0e-14)

    A25circ = HRR.rr_harmonic_accel(Xcirc, Vcirc, q; include25=true, include35=false)
    @test dot(A25circ, Vcirc) < 0.0
    @test isapprox(dot(A25circ, Xcirc), 0.0; atol=1.0e-14, rtol=1.0e-10)

    # The RR acceleration is linear in q at leading EMRI order, but the
    # coefficients still contain u = m / R = (1 + q) / R, so a finite-q test
    # should only be approximately linear. Using smaller q values keeps this a
    # genuine small-mass-ratio check without over-constraining higher-order
    # q-dependence.
    Aq1 = HRR.rr_harmonic_accel(X, V, 1.0e-7; include25=true, include35=true)
    Aq2 = HRR.rr_harmonic_accel(X, V, 2.0e-7; include25=true, include35=true)
    @test isapprox(norm(Aq2) / norm(Aq1), 2.0; atol=5.0e-6)

    @test_throws DomainError HRR.rr_harmonic_accel(zeros(3), V, q)
end

@testset "BT to harmonic correction formulas" begin
    X = @SVector [8.1, -3.2, 4.7]
    V = @SVector [0.07, -0.03, 0.02]
    q = 1.0e-5

    m = 1.0 * (1 + q)
    ν = q / (1 + q)^2
    R = norm(X)
    N = X / R
    Rdot = dot(N, V)
    V2 = dot(V, V)
    u = m / R

    A25H = HRR.rr_harmonic_accel(X, V, q; include25=true, include35=false)
    A25BT = HRR.rr_bt_accel(X, V, q; include25=true, include35=false)
    Δ25 = HRR.rr_bt_to_harmonic_correction(X, V, q; include25=true, include35=false)

    pref25BT = (8 / 5) * m^2 * ν / R^3
    A5 = 18 * V2 + (2 / 3) * u - 25 * Rdot^2
    B5 = 6 * V2 - 2 * u - 15 * Rdot^2
    expected25BT = pref25BT * (Rdot * A5 * N - B5 * V)

    @test A25BT ≈ expected25BT

    expectedΔ25 = pref25BT * (
        ((-15 * V2 + 5 * u + 25 * Rdot^2) * Rdot) * N -
        (-5 * V2 + 5 * u + 15 * Rdot^2) * V
    )

    @test Δ25 ≈ expectedΔ25
    @test Δ25 ≈ A25H - A25BT
    @test_throws DomainError HRR.rr_bt_accel(zeros(3), V, q)

    A35BT = HRR.rr_bt_accel(X, V, q; include25=false, include35=true)
    pref35BT = (8 / 5) * m^2 * ν / R^3
    A7 =
        (87 / 14 - 48 * ν) * V2^2 -
        (5379 / 28 + (136 / 3) * ν) * V2 * u +
        (25 / 2) * (1 + 5 * ν) * V2 * Rdot^2 +
        (1353 / 4 + 133 * ν) * Rdot^2 * u -
        (35 / 2) * (1 - ν) * Rdot^4 +
        (160 / 7 + (55 / 3) * ν) * u^2
    B7 =
        -(27 / 14) * V2^2 -
        (4861 / 84 + (58 / 3) * ν) * V2 * u +
        (3 / 2) * (13 - 37 * ν) * V2 * Rdot^2 +
        (2591 / 12 + 97 * ν) * Rdot^2 * u -
        (25 / 2) * (1 - 7 * ν) * Rdot^4 +
        (1 / 3) * (776 / 7 + 55 * ν) * u^2
    expected35BT = pref35BT * (Rdot * A7 * N - B7 * V)

    @test A35BT ≈ expected35BT
end

@testset "Harmonic RR lift and transform" begin
    Acoord = @SVector [0.01, -0.02, 0.03]
    Vflat = @SVector [0.1, -0.05, 0.07]
    gflat = [-1.0 0.0 0.0 0.0;
             0.0 1.0 0.0 0.0;
             0.0 0.0 1.0 0.0;
             0.0 0.0 0.0 1.0]

    Γflat = HRR.gamma_from_metric_velocity(gflat, Vflat)
    a4_flat = HRR.lift_coord_accel_to_four(Acoord, Vflat, gflat)
    uflat = @SVector [Γflat, Γflat * Vflat[1], Γflat * Vflat[2], Γflat * Vflat[3]]
    uflat_cov = gflat * uflat

    @test isapprox(dot(uflat_cov, a4_flat), 0.0; atol=1.0e-12)
    @test HRR.coord_accel_from_four(a4_flat, Vflat, Γflat) ≈ Acoord

    spin = 0.3
    q = 1.0e-5
    xBL = @SVector [25.0, π / 3, π / 4]
    xH = iChimera.HarmonicCoords.xBLtoH(xBL, spin)
    Vweak = @SVector [1.0e-3, -7.0e-4, 5.0e-4]
    gH = iChimera.HarmonicCoords.g_μν_H(xH, spin)

    Γweak = HRR.gamma_from_metric_velocity(gH, Vweak)
    a4_weak = HRR.lift_coord_accel_to_four(Acoord, Vweak, gH)
    uweak = @SVector [Γweak, Γweak * Vweak[1], Γweak * Vweak[2], Γweak * Vweak[3]]
    uweak_cov = gH * uweak

    @test isapprox(dot(uweak_cov, a4_weak), 0.0; atol=1.0e-11, rtol=1.0e-11)
    @test HRR.coord_accel_from_four(a4_weak, Vweak, Γweak) ≈ Acoord

    trial_a4_H = @SVector [0.12, 0.03, -0.02, 0.01]
    trial_a4_BL = HRR.four_accel_H_to_BL(trial_a4_H, xH, spin)
    @test trial_a4_BL[1] == trial_a4_H[1]
    @test trial_a4_BL[2:4] ≈ iChimera.HarmonicCoords.aHtoBL(xH, zeros(3), collect(trial_a4_H[2:4]), spin)
    @test HRR.four_accel_coord_to_BL(trial_a4_H, xH, spin; coord_system=:harmonic) ≈ trial_a4_BL

    corr_BL = HRR.rr_bt_to_harmonic_correction_four_accel_BL(xH, Vweak, spin, q; include25=true, include35=true)
    Acoord_corr = HRR.rr_bt_to_harmonic_correction(xH, Vweak, q; include25=true, include35=true)
    a4_corr_H = HRR.lift_coord_accel_to_four(Acoord_corr, Vweak, gH)
    @test corr_BL ≈ HRR.four_accel_H_to_BL(a4_corr_H, xH, spin)

    chimera_BL = @SVector [0.01, -0.02, 0.03, -0.04]
    fixed_BL = MVector{4, Float64}(undef)
    scratch_corr_BL = MVector{4, Float64}(undef)
    scratch_corr_H = MVector{4, Float64}(undef)
    scratch_Acoord_corr = MVector{3, Float64}(undef)
    scratch_Acoord_BT = MVector{3, Float64}(undef)
    HRR.add_bt_to_harmonic_correction_four_accel_BL!(
        fixed_BL,
        chimera_BL,
        scratch_corr_BL,
        scratch_corr_H,
        scratch_Acoord_corr,
        scratch_Acoord_BT,
        xH,
        Vweak,
        spin,
        q,
    )
    @test fixed_BL ≈ chimera_BL + corr_BL

    xC = @SVector [25.0, -6.0, 8.0]
    Vcart = @SVector [9.0e-4, -6.0e-4, 4.0e-4]
    gC = iChimera.CartesianCoords.g_μν_H(xC, spin)

    xC_from_BL = iChimera.CartesianCoords.xBLtoH(xBL, spin)
    xC_from_BL_mut = MVector{3, Float64}(undef)
    iChimera.CartesianCoords.xBLtoH!(xC_from_BL_mut, collect(xBL), spin)
    @test xC_from_BL ≈ xC_from_BL_mut

    Γcart = HRR.gamma_from_metric_velocity(gC, Vcart)
    a4_cart = HRR.lift_coord_accel_to_four(Acoord, Vcart, gC)
    ucart = @SVector [Γcart, Γcart * Vcart[1], Γcart * Vcart[2], Γcart * Vcart[3]]
    ucart_cov = gC * ucart

    @test isapprox(dot(ucart_cov, a4_cart), 0.0; atol=1.0e-11, rtol=1.0e-11)
    @test HRR.coord_accel_from_four(a4_cart, Vcart, Γcart) ≈ Acoord

    trial_a4_C = @SVector [0.08, -0.01, 0.04, -0.03]
    trial_a4_BL_cart = HRR.four_accel_coord_to_BL(trial_a4_C, xC, spin; coord_system=:cartesian)
    @test trial_a4_BL_cart[1] == trial_a4_C[1]
    @test trial_a4_BL_cart[2:4] ≈ iChimera.CartesianCoords.aHtoBL(xC, zeros(3), collect(trial_a4_C[2:4]), spin)
end

@testset "RR model plumbing" begin
    @test iChimera.Inspiral.normalize_rr_model(:chimera) === :chimera
    @test iChimera.Inspiral.normalize_rr_model("harmonic_pn") === :harmonic_pn
    @test iChimera.Inspiral.normalize_rr_model(:chimera_bt_to_harmonic) === :chimera_bt_to_harmonic
    @test iChimera.Inspiral.normalize_rr_model("chimera_bt_to_harmonic") === :chimera_bt_to_harmonic
    @test_throws ArgumentError iChimera.Inspiral.normalize_rr_model(:not_a_model)

    chimera_fname = iChimera.Inspiral.solution_fname(0.3, 12.0, 0.2, π / 3, 1.0e-5, 0.1, 0.2, 0.3, 3, 2, "harmonic", "/tmp/"; rr_model=:chimera)
    harmonic_fname = iChimera.Inspiral.solution_fname(0.3, 12.0, 0.2, π / 3, 1.0e-5, 0.1, 0.2, 0.3, 3, 2, "harmonic", "/tmp/"; rr_model=:harmonic_pn)
    corrected_fname = iChimera.Inspiral.solution_fname(0.3, 12.0, 0.2, π / 3, 1.0e-5, 0.1, 0.2, 0.3, 3, 2, "harmonic", "/tmp/"; rr_model=:chimera_bt_to_harmonic)
    harmonic_wave = iChimera.Inspiral.waveform_fname(0.3, 12.0, 0.2, π / 3, 1.0e-5, 0.1, 0.2, 0.3, 1.0, 20.0, 10.0, 3, 2, 4, 3, "harmonic", "/tmp/"; rr_model=:harmonic_pn)
    corrected_wave = iChimera.Inspiral.waveform_fname(0.3, 12.0, 0.2, π / 3, 1.0e-5, 0.1, 0.2, 0.3, 1.0, 20.0, 10.0, 3, 2, 4, 3, "harmonic", "/tmp/"; rr_model=:chimera_bt_to_harmonic)

    @test !occursin("_rrharmonicpn", chimera_fname)
    @test occursin("_rrharmonicpn.h5", harmonic_fname)
    @test occursin("_rrchimeraharmonicfix.h5", corrected_fname)
    @test occursin("_rrharmonicpn.h5", harmonic_wave)
    @test occursin("_rrchimeraharmonicfix.h5", corrected_wave)
    @test chimera_fname != harmonic_fname
    @test corrected_fname != chimera_fname
    @test corrected_fname != harmonic_fname

    @test iChimera.Inspiral.coordinate_module("harmonic") === iChimera.HarmonicCoords
    @test iChimera.Inspiral.coordinate_module("cartesian") === iChimera.CartesianCoords
    @test_throws ArgumentError iChimera.Inspiral.coordinate_module("not_a_coord")
end
