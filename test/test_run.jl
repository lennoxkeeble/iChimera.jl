include("test_params.jl");


@time iChimera.compute_inspiral(emri);
@time iChimera.compute_waveform(emri);

# load trajectory
t, r, θ, ϕ = iChimera.load_trajectory(emri);
t = t * MtoSecs;

# load fluxes
t_Fluxes, EE, LL, QQ, CC, pArray, ecc, θminArray = iChimera.load_constants_of_motion(emri);
t_Fluxes, Edot, Ldot, Qdot, Cdot = iChimera.load_fluxes(emri);
t_Fluxes = t_Fluxes * MtoSecs;

# compute iota
ι = @. acos(LL / sqrt(LL^2 + CC));

t_wf, h_plus, h_cross = iChimera.load_waveform(emri);
h_plus, h_cross = h_plus * strain_to_SI, h_cross * strain_to_SI;
t_wf *= MtoSecs;