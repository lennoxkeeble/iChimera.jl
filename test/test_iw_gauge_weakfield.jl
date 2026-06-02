include(joinpath(@__DIR__, "iw_gauge_diagnostic_helpers.jl"))

function weakfield_case(R::Float64, q::Float64, spin::Float64, eps_fd::Float64, gauge_sign::Float64; zero_potentials::Bool)
    xH = [R, 0.0, 0.0]
    vH = [0.006, 0.025, 0.004]
    aH = -xH / norm(xH)^3
    jH = newtonian_jerk(xH, vH)
    xBL = iChimera.HarmonicCoords.xHtoBL(xH, spin)
    rH = norm(xH)
    v = SH.norm_3d(vH)

    rr = zero_potentials ? zero_rr_data() : dummy_rr_data(scale=1.0e-1)
    aSF_H_bt = zeros(4)
    aSF_BL_bt = zeros(4)
    aSF_H_target = zeros(4)
    aSF_BL_target = zeros(4)

    SHIW.aRRα(
        aSF_H_bt,
        aSF_BL_bt,
        xH,
        v,
        vH,
        aH,
        jH,
        xBL,
        rH,
        spin,
        q,
        rr.Vrr,
        rr.dVrr_dt,
        rr.Virr,
        rr.dVrr_da,
        rr.dVirr_dt,
        rr.dVirr_da;
        alpha_current=4.0,
        beta_current=5.0,
        alpha_target=4.0,
        beta_target=5.0,
        gauge_sign=gauge_sign,
        eps_fd=eps_fd,
    )

    SHIW.aRRα(
        aSF_H_target,
        aSF_BL_target,
        xH,
        v,
        vH,
        aH,
        jH,
        xBL,
        rH,
        spin,
        q,
        rr.Vrr,
        rr.dVrr_dt,
        rr.Virr,
        rr.dVrr_da,
        rr.dVirr_dt,
        rr.dVirr_da;
        alpha_current=4.0,
        beta_current=5.0,
        alpha_target=0.0,
        beta_target=0.0,
        gauge_sign=gauge_sign,
        eps_fd=eps_fd,
    )

    δa_coord = self_force_difference_to_coord_accel(aSF_H_target, aSF_H_bt, vH, xH, spin)
    expected = expected_IW_shift(xH, vH, q, 4.0, 5.0, 0.0, 0.0)
    rel_err = norm(δa_coord - expected) / max(norm(expected), eps(Float64))

    return (
        R=R,
        q=q,
        eps_fd=eps_fd,
        gauge_sign=gauge_sign,
        zero_potentials=zero_potentials,
        δa_coord=δa_coord,
        expected=expected,
        rel_err=rel_err,
        abs_err=norm(δa_coord - expected),
    )
end

function run_iw_gauge_weakfield_diagnostic(;
    spin::Float64=1.0e-3,
    q_values::Vector{Float64}=[1.0e-4, 1.0e-5, 1.0e-6],
    radii::Vector{Float64}=[300.0, 600.0, 1200.0],
    eps_scales::Vector{Float64}=[1.0e-4, 1.0e-5, 1.0e-6],
)
    println("Weak-field Iyer-Will benchmark:")
    println("  spin = $(spin)")
    println("  target map = (4,5) -> (0,0)")
    println()

    sign_results = Dict{Float64, Vector{NamedTuple}}()
    for gauge_sign in (-1.0, 1.0)
        cases = NamedTuple[]
        for R in radii
            eps_fd = 1.0e-5 * max(1.0, R)
            push!(cases, weakfield_case(R, 1.0e-5, spin, eps_fd, gauge_sign; zero_potentials=false))
        end
        sign_results[gauge_sign] = cases
    end

    mean_err(sign) = mean(case.rel_err for case in sign_results[sign])
    best_sign = mean_err(-1.0) <= mean_err(1.0) ? -1.0 : 1.0

    println("CHECK 6 sign comparison:")
    @printf("  gauge_sign = -1 mean relative mismatch: %.3e\n", mean_err(-1.0))
    @printf("  gauge_sign = +1 mean relative mismatch: %.3e\n", mean_err(1.0))
    @printf("  preferred gauge_sign from weak-field benchmark: %+0.1f\n", best_sign)
    println()

    chosen_cases = sign_results[best_sign]
    println("CHECK 6 radius scaling with dummy nonzero RR potentials:")
    for case in chosen_cases
        @printf("  R = %7.1f, rel mismatch = %.3e, |delta a - expected| = %.3e\n", case.R, case.rel_err, case.abs_err)
    end
    println()

    println("CHECK 7 pure-gauge isolated test (all RR potentials zero):")
    pure_gauge_cases = [
        weakfield_case(R, 1.0e-5, spin, 1.0e-5 * max(1.0, R), best_sign; zero_potentials=true)
        for R in radii
    ]
    for case in pure_gauge_cases
        @printf("  R = %7.1f, rel mismatch = %.3e, |delta a - expected| = %.3e\n", case.R, case.rel_err, case.abs_err)
    end
    println()

    println("CHECK 9 finite-difference stability at q = 1e-5:")
    fd_cases = [
        weakfield_case(radii[end], 1.0e-5, spin, eps_scale * max(1.0, radii[end]), best_sign; zero_potentials=true)
        for eps_scale in eps_scales
    ]
    for case in fd_cases
        @printf("  eps_fd/R = %.1e, rel mismatch = %.3e, |delta a| = %.3e\n", case.eps_fd / case.R, case.rel_err, norm(case.δa_coord))
    end
    println()

    println("CHECK 13 q-scaling at R = $(radii[end]):")
    q_cases = [
        weakfield_case(radii[end], q, spin, 1.0e-5 * max(1.0, radii[end]), best_sign; zero_potentials=true)
        for q in q_values
    ]
    ref_q = q_cases[2]
    ref_norm = norm(ref_q.δa_coord)
    ref_q_value = ref_q.q
    for case in q_cases
        ratio_measured = norm(case.δa_coord) / ref_norm
        ratio_expected = case.q / ref_q_value
        @printf("  q = %.1e, measured ratio = %.3e, expected ratio = %.3e\n", case.q, ratio_measured, ratio_expected)
    end

    return (
        best_sign=best_sign,
        sign_results=sign_results,
        pure_gauge_cases=pure_gauge_cases,
        fd_cases=fd_cases,
        q_cases=q_cases,
    )
end

results_weakfield = run_iw_gauge_weakfield_diagnostic()

