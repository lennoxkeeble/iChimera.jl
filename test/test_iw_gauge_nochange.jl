include(joinpath(@__DIR__, "iw_gauge_diagnostic_helpers.jl"))

function run_iw_gauge_nochange_diagnostic(;
    spin::Float64=1.0e-2,
    p::Float64=14.0,
    e::Float64=0.35,
    θmin::Float64=1.1,
    sign_Lz::Int=1,
    q::Float64=1.0e-5,
    sample_fraction::Float64=0.37,
)
    print_repo_summary()

    setup = build_geodesic_integrator(
        a=spin,
        p=p,
        e=e,
        θmin=θmin,
        sign_Lz=sign_Lz,
        nOrbits=1.0,
    )
    target_t = sample_fraction * setup.t_geodesic
    step_integrator_to!(setup.integrator, target_t)

    ws = allocate_flux_workspace()
    local_state, old_flux, new_flux = compute_local_flux_pair!(
        ws,
        setup.integrator,
        spin,
        p,
        e,
        θmin,
        q,
        setup.E,
        setup.L,
        setup.Q,
        setup.C,
        setup.p3,
        setup.p4,
        setup.zp,
        setup.zm;
        onepn=0.0,
        twopn=0.0,
        twop5pn=0.0,
        v2power=0,
        alpha_current=4.0,
        beta_current=5.0,
        alpha_target=4.0,
        beta_target=5.0,
        gauge_sign=-1.0,
    )

    abs_H, rel_H = abs_rel_diff(ws.aSF_H_old, ws.aSF_H_new)
    abs_BL, rel_BL = abs_rel_diff(ws.aSF_BL_old, ws.aSF_BL_new)
    abs_flux_E, rel_flux_E = abs_rel_diff([old_flux.Edot], [new_flux.Edot])
    abs_flux_L, rel_flux_L = abs_rel_diff([old_flux.Ldot], [new_flux.Ldot])
    abs_flux_Q, rel_flux_Q = abs_rel_diff([old_flux.Qdot], [new_flux.Qdot])

    xH_local = [200.0, 30.0, 10.0]
    vH_local = [0.015, 0.055, 0.011]
    aH_local = -xH_local / norm(xH_local)^3
    jH_local = newtonian_jerk(xH_local, vH_local)

    ξ, dtξ, dxξ, dttξ, dtdxξ = SHIW._gauge_data(
        xH_local,
        vH_local,
        aH_local,
        jH_local,
        q;
        alpha_current=4.0,
        beta_current=5.0,
        alpha_target=4.0,
        beta_target=5.0,
        gauge_sign=-1.0,
    )

    Lξg = SHIW._lie_xi_gK(
        xH_local,
        vH_local,
        aH_local,
        jH_local,
        q,
        spin;
        alpha_current=4.0,
        beta_current=5.0,
        alpha_target=4.0,
        beta_target=5.0,
        gauge_sign=-1.0,
    )

    rr_dummy = dummy_rr_data()
    h0 = SHIW._base_hRR(rr_dummy.Vrr, rr_dummy.Virr)
    dh0 = SHIW._base_dhRR(rr_dummy.dVrr_dt, rr_dummy.dVrr_da, rr_dummy.dVirr_dt, rr_dummy.dVirr_da)
    h_tilde = SHIW._hRR_gauged(
        xH_local,
        vH_local,
        aH_local,
        jH_local,
        q,
        spin,
        rr_dummy.Vrr,
        rr_dummy.Virr;
        alpha_current=4.0,
        beta_current=5.0,
        alpha_target=4.0,
        beta_target=5.0,
        gauge_sign=-1.0,
    )
    dh_tilde = SHIW._dhRR_gauged(
        xH_local,
        vH_local,
        aH_local,
        jH_local,
        q,
        spin,
        rr_dummy.dVrr_dt,
        rr_dummy.dVrr_da,
        rr_dummy.dVirr_dt,
        rr_dummy.dVirr_da;
        alpha_current=4.0,
        beta_current=5.0,
        alpha_target=4.0,
        beta_target=5.0,
        gauge_sign=-1.0,
    )

    results = (
        sample_time=target_t,
        Ω=setup.Ω,
        T_max=setup.T_max,
        abs_H=abs_H,
        rel_H=rel_H,
        abs_BL=abs_BL,
        rel_BL=rel_BL,
        abs_flux_E=abs_flux_E,
        rel_flux_E=rel_flux_E,
        abs_flux_L=abs_flux_L,
        rel_flux_L=rel_flux_L,
        abs_flux_Q=abs_flux_Q,
        rel_flux_Q=rel_flux_Q,
        norm_xi=norm(ξ),
        norm_dtxi=norm(dtξ),
        norm_dxxi=norm(dxξ),
        norm_dttxi=norm(dttξ),
        norm_dtdxxi=norm(dtdxξ),
        norm_Lxig=norm(Lξg),
        h_match_norm=norm(h_tilde - h0),
        dh_match_norm=norm(dh_tilde - dh0),
    )

    println()
    @printf("CHECK 2 no-gauge self-force H: abs = %.3e, rel = %.3e\n", results.abs_H, results.rel_H)
    @printf("CHECK 2 no-gauge self-force BL: abs = %.3e, rel = %.3e\n", results.abs_BL, results.rel_BL)
    @printf("CHECK 2 no-gauge flux dE/dt: abs = %.3e, rel = %.3e\n", results.abs_flux_E, results.rel_flux_E)
    @printf("CHECK 2 no-gauge flux dL/dt: abs = %.3e, rel = %.3e\n", results.abs_flux_L, results.rel_flux_L)
    @printf("CHECK 2 no-gauge flux dQ/dt: abs = %.3e, rel = %.3e\n", results.abs_flux_Q, results.rel_flux_Q)
    @printf("CHECK 3 gauge-vector norms: |xi| = %.3e, |dt xi| = %.3e, |dx xi| = %.3e\n", results.norm_xi, results.norm_dtxi, results.norm_dxxi)
    @printf("CHECK 3 higher derivatives: |dtt xi| = %.3e, |dt dx xi| = %.3e\n", results.norm_dttxi, results.norm_dtdxxi)
    @printf("CHECK 4 Lie derivative norm: %.3e\n", results.norm_Lxig)
    @printf("CHECK 5 h equality norms: |h_tilde - h| = %.3e, |dh_tilde - dh| = %.3e\n", results.h_match_norm, results.dh_match_norm)

    return results
end

results_nochange = run_iw_gauge_nochange_diagnostic()

