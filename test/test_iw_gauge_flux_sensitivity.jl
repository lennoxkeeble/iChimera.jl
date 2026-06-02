include(joinpath(@__DIR__, "iw_gauge_diagnostic_helpers.jl"))

function run_iw_flux_sensitivity(;
    spin::Float64=0.1,
    p::Float64=11.5,
    e::Float64=0.45,
    θmin::Float64=1.05,
    sign_Lz::Int=1,
    q::Float64=1.0e-5,
    nOrbits::Float64=2.0,
    n_samples::Int=200,
    alpha_current::Float64=4.0,
    beta_current::Float64=5.0,
    alpha_target::Float64=0.0,
    beta_target::Float64=0.0,
    gauge_sign::Float64=-1.0,
)
    setup = build_geodesic_integrator(
        a=spin,
        p=p,
        e=e,
        θmin=θmin,
        sign_Lz=sign_Lz,
        nOrbits=nOrbits,
    )
    sample_times = collect(range(0.0, setup.t_geodesic, length=n_samples))
    ws = allocate_flux_workspace()

    Edot_old = zeros(n_samples)
    Ldot_old = zeros(n_samples)
    Qdot_old = zeros(n_samples)
    Edot_new = zeros(n_samples)
    Ldot_new = zeros(n_samples)
    Qdot_new = zeros(n_samples)

    for (idx, t_sample) in enumerate(sample_times)
        step_integrator_to!(setup.integrator, t_sample)
        _, old_flux, new_flux = compute_local_flux_pair!(
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
            alpha_current=alpha_current,
            beta_current=beta_current,
            alpha_target=alpha_target,
            beta_target=beta_target,
            gauge_sign=gauge_sign,
        )
        Edot_old[idx] = old_flux.Edot
        Ldot_old[idx] = old_flux.Ldot
        Qdot_old[idx] = old_flux.Qdot
        Edot_new[idx] = new_flux.Edot
        Ldot_new[idx] = new_flux.Ldot
        Qdot_new[idx] = new_flux.Qdot
    end

    metrics = (
        Edot_max_abs=maximum(abs.(Edot_new .- Edot_old)),
        Ldot_max_abs=maximum(abs.(Ldot_new .- Ldot_old)),
        Qdot_max_abs=maximum(abs.(Qdot_new .- Qdot_old)),
        Edot_mean_old=mean(Edot_old),
        Ldot_mean_old=mean(Ldot_old),
        Qdot_mean_old=mean(Qdot_old),
        Edot_mean_new=mean(Edot_new),
        Ldot_mean_new=mean(Ldot_new),
        Qdot_mean_new=mean(Qdot_new),
    )

    println("Fixed-geodesic flux sensitivity test:")
    @printf("  orbit: spin=%.3f, p=%.3f, e=%.3f, θmin=%.3f, nOrbits=%.1f\n", spin, p, e, θmin, nOrbits)
    @printf("  gauge map: (%.1f, %.1f) -> (%.1f, %.1f), gauge_sign=%+.1f\n", alpha_current, beta_current, alpha_target, beta_target, gauge_sign)
    println("  BL-time frequencies = $(setup.Ω)")
    @printf("  T_max = %.6f, duration = %.6f\n", setup.T_max, setup.t_geodesic)
    @printf("  max |ΔEdot| = %.3e\n", metrics.Edot_max_abs)
    @printf("  max |ΔLdot| = %.3e\n", metrics.Ldot_max_abs)
    @printf("  max |ΔQdot| = %.3e\n", metrics.Qdot_max_abs)
    @printf("  <Edot> old/new = %.12e / %.12e\n", metrics.Edot_mean_old, metrics.Edot_mean_new)
    @printf("  <Ldot> old/new = %.12e / %.12e\n", metrics.Ldot_mean_old, metrics.Ldot_mean_new)
    @printf("  <Qdot> old/new = %.12e / %.12e\n", metrics.Qdot_mean_old, metrics.Qdot_mean_new)

    return (
        times=sample_times,
        Edot_old=Edot_old,
        Ldot_old=Ldot_old,
        Qdot_old=Qdot_old,
        Edot_new=Edot_new,
        Ldot_new=Ldot_new,
        Qdot_new=Qdot_new,
        metrics=metrics,
        Ω=setup.Ω,
        T_max=setup.T_max,
        t_geodesic=setup.t_geodesic,
    )
end

results_flux_sensitivity = run_iw_flux_sensitivity()
