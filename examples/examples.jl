include("params.jl");
@time iChimera.compute_inspiral(emri; JIT = true);
@time iChimera.compute_inspiral(emri; JIT = false);
@time iChimera.compute_waveform(emri);