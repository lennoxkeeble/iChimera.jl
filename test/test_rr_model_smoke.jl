using iChimera
using Test
using HDF5

# Manual smoke tests for the inspiral pipeline.
# These are intentionally separate from the unit-test file because they execute
# short inspirals and therefore cost much more than the local formula tests.

function smoke_emri(path::String; a::Float64, coordinates::String)
    p = 10.0
    e = 0.2
    inclination = 30.0
    inclination_type = "iota"
    sign_Lz = 1
    mass_ratio = 1.0e-5

    lmax_mass_fluxes = 3
    lmax_current_fluxes = 2
    lmax_mass_waveform = 4
    lmax_current_waveform = 3

    OnePN = 1.0
    TwoPN = 1.0
    TwoPointFivePN = 1.0

    psi0 = 0.1
    chi0 = 0.2
    phi0 = 0.3

    frame = "Source"
    ThetaS = 10.0
    PhiS = 20.0
    ThetaK = 30.0
    PhiK = 40.0
    ThetaObs = 50.0
    PhiObs = 60.0

    dt_save = 5.0
    T_secs = 50.0
    M = 1.0e6 * iChimera.Msol
    reltol = 1.0e-12
    abstol = 1.0e-12
    compute_SF_frac = 0.05
    save_every = 8

    save_traj = false
    save_constants = true
    save_fluxes = true
    save_gamma = false

    return iChimera.EMRI(
        a, p, e, inclination, inclination_type, sign_Lz, mass_ratio,
        lmax_mass_fluxes, lmax_current_fluxes,
        lmax_mass_waveform, lmax_current_waveform,
        coordinates, OnePN, TwoPN, TwoPointFivePN,
        psi0, chi0, phi0,
        frame, ThetaS, PhiS, ThetaK, PhiK, ThetaObs, PhiObs,
        dt_save, path, T_secs, M, reltol, abstol, compute_SF_frac,
        save_every, save_traj, save_constants, save_fluxes, save_gamma
    )
end

@testset "Manual rr_model smoke tests" begin
    mktempdir() do tmpdir
        chimera_path = joinpath(tmpdir, "chimera") * "/"
        harmonic_path = joinpath(tmpdir, "harmonic_pn") * "/"
        harmonic_cart_path = joinpath(tmpdir, "harmonic_pn_cartesian") * "/"
        corrected_path = joinpath(tmpdir, "chimera_bt_to_harmonic") * "/"
        mkpath(chimera_path)
        mkpath(harmonic_path)
        mkpath(harmonic_cart_path)
        mkpath(corrected_path)

        chimera_emri = smoke_emri(chimera_path; a=0.3, coordinates="cartesian")
        harmonic_emri = smoke_emri(harmonic_path; a=0.7, coordinates="harmonic")
        harmonic_cart_emri = smoke_emri(harmonic_cart_path; a=0.7, coordinates="cartesian")
        corrected_emri = smoke_emri(corrected_path; a=0.7, coordinates="harmonic")

        iChimera.compute_inspiral(chimera_emri; rr_model=:chimera)
        iChimera.compute_inspiral(harmonic_emri; rr_model=:harmonic_pn)
        iChimera.compute_inspiral(harmonic_cart_emri; rr_model=:harmonic_pn)
        iChimera.compute_inspiral(corrected_emri; rr_model=:chimera_bt_to_harmonic)

        chimera_sol = iChimera.Inspiral.solution_fname(
            chimera_emri.a, chimera_emri.p, chimera_emri.e, chimera_emri.θmin,
            chimera_emri.mass_ratio, chimera_emri.psi0, chimera_emri.chi0, chimera_emri.phi0,
            chimera_emri.lmax_mass_fluxes, chimera_emri.lmax_current_fluxes,
            chimera_emri.coordinates, chimera_emri.path; rr_model=:chimera
        )
        harmonic_sol = iChimera.Inspiral.solution_fname(
            harmonic_emri.a, harmonic_emri.p, harmonic_emri.e, harmonic_emri.θmin,
            harmonic_emri.mass_ratio, harmonic_emri.psi0, harmonic_emri.chi0, harmonic_emri.phi0,
            harmonic_emri.lmax_mass_fluxes, harmonic_emri.lmax_current_fluxes,
            harmonic_emri.coordinates, harmonic_emri.path; rr_model=:harmonic_pn
        )
        harmonic_cart_sol = iChimera.Inspiral.solution_fname(
            harmonic_cart_emri.a, harmonic_cart_emri.p, harmonic_cart_emri.e, harmonic_cart_emri.θmin,
            harmonic_cart_emri.mass_ratio, harmonic_cart_emri.psi0, harmonic_cart_emri.chi0, harmonic_cart_emri.phi0,
            harmonic_cart_emri.lmax_mass_fluxes, harmonic_cart_emri.lmax_current_fluxes,
            harmonic_cart_emri.coordinates, harmonic_cart_emri.path; rr_model=:harmonic_pn
        )
        corrected_sol = iChimera.Inspiral.solution_fname(
            corrected_emri.a, corrected_emri.p, corrected_emri.e, corrected_emri.θmin,
            corrected_emri.mass_ratio, corrected_emri.psi0, corrected_emri.chi0, corrected_emri.phi0,
            corrected_emri.lmax_mass_fluxes, corrected_emri.lmax_current_fluxes,
            corrected_emri.coordinates, corrected_emri.path; rr_model=:chimera_bt_to_harmonic
        )

        @test isfile(chimera_sol)
        @test isfile(harmonic_sol)
        @test isfile(harmonic_cart_sol)
        @test isfile(corrected_sol)
        @test chimera_sol != harmonic_sol
        @test harmonic_cart_sol != harmonic_sol
        @test chimera_sol != corrected_sol
        @test harmonic_sol != corrected_sol
        @test occursin("_rrharmonicpn", harmonic_sol)
        @test occursin("_rrharmonicpn", harmonic_cart_sol)
        @test occursin("_rrchimeraharmonicfix", corrected_sol)

        _, Edot_chimera, Ldot_chimera, Qdot_chimera, Cdot_chimera = iChimera.load_fluxes(chimera_emri; rr_model=:chimera)
        _, Edot_harmonic, Ldot_harmonic, Qdot_harmonic, Cdot_harmonic = iChimera.load_fluxes(harmonic_emri; rr_model=:harmonic_pn)
        _, Edot_harmonic_cart, Ldot_harmonic_cart, Qdot_harmonic_cart, Cdot_harmonic_cart = iChimera.load_fluxes(harmonic_cart_emri; rr_model=:harmonic_pn)
        _, Edot_corrected, Ldot_corrected, Qdot_corrected, Cdot_corrected = iChimera.load_fluxes(corrected_emri; rr_model=:chimera_bt_to_harmonic)

        @test all(isfinite, Edot_chimera)
        @test all(isfinite, Ldot_chimera)
        @test all(isfinite, Qdot_chimera)
        @test all(isfinite, Cdot_chimera)

        @test all(isfinite, Edot_harmonic)
        @test all(isfinite, Ldot_harmonic)
        @test all(isfinite, Qdot_harmonic)
        @test all(isfinite, Cdot_harmonic)

        @test all(isfinite, Edot_harmonic_cart)
        @test all(isfinite, Ldot_harmonic_cart)
        @test all(isfinite, Qdot_harmonic_cart)
        @test all(isfinite, Cdot_harmonic_cart)

        @test all(isfinite, Edot_corrected)
        @test all(isfinite, Ldot_corrected)
        @test all(isfinite, Qdot_corrected)
        @test all(isfinite, Cdot_corrected)
    end
end
