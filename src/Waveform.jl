#=

    In this module we project the metric perturbation from the kludge scheme in arXiv:1109.0572v2 into the TT gauge.

=#

module Waveform
using Combinatorics, LinearAlgebra

"""
# Common Arguments in this module
- `r::Float64`: observer distance.
- `Θ::Float64`: observer polar orientation.
- `ϕ::Float64`: observer azimuthal orientation.
- `Mij2::AbstractArray`: second derivative of the mass quadrupole (Eq. 48).
- `Mijk3::AbstractArray`: third derivative of the mass quadrupole (Eq. 48).
- `Mijkl4::AbstractArray`: fourth derivative of the mass quadrupole (Eq. 85).
- `Sij2::AbstractArray`: second derivative of the current quadrupole (Eq. 49).
- `Sijk3::AbstractArray`: third derivative of the current quadrupole (Eq. 86).
- `ThetaObs::Float64`: polar angle of the observer in the source frame (i.e., Cartesian coordinate system with origin at the center of the BH)
- `PhiObs::Float64`: azimuthal angle of the observer in the source frame (i.e., Cartesian coordinate system with origin at the center of the BH)
- `Phi0::Float64`: azimuthal angle of the underlying geodesic orbit at t=0 (assuming that in the specified hplus, hcross, the orbit was started instead with phi(t=0) = 0 automatically)
- `ThetaSource::Float64`: polar angle of the unit vector pointing from the Solar System Barycenter (SSB) to the EMRI system in the SSB frame
- `PhiSource::Float64`: azimuthal angle of the unit vector pointing to the EMRI system in the SSB frame
- `ThetaKerr::Float64`: polar angle of the MBH's spin vector in the SSB frame
- `PhiKerr::Float64`: azimuthal angle of the MBH's spin vector in the SSB frame
"""

const spatial_indices_3::Array = [[x, y, z] for x=1:3, y=1:3, z=1:3]
const εkl::Array{Vector} = [[levicivita(spatial_indices_3[k, l, i]) for i = 1:3] for k=1:3, l=1:3]

δ(i::Int, j::Int) = i == j ? 1.0 : 0.0
@inline outer(x::Vector{Float64}, y::Vector{Float64}) = [x[i] * y[j] for i in eachindex(x), j in eachindex(y)]


# returns plus and cross polarized waveforms with respect to the polarization vectors x_w, y_w as defined in Eq. 6 in https://arxiv.org/pdf/2104.04582
function compute_wave_polarizations(nPoints::Int, r::Float64, ThetaSource::Float64, PhiSource::Float64, ThetaKerr::Float64, PhiKerr::Float64, Mij2::AbstractArray, Mijk3::AbstractArray, Mijkl4::AbstractArray, Sij2::AbstractArray, Sijk3::AbstractArray, mass_ratio::Float64)
    hij = [zeros(nPoints) for i=1:3, j=1:3];
    hij_TT = [zeros(nPoints) for i=1:3, j=1:3];
    hplus = zeros(nPoints);
    hcross = zeros(nPoints);
    
    # R_ssb ≡ unit vector pointing from solar system barycenter (SSB) in the direction of the EMRI system in the SSB frame
    R_ssb = [sin(ThetaSource) * cos(PhiSource), sin(ThetaSource) * sin(PhiSource), cos(ThetaSource)]
    
    # S_ssb ≡ unit vector pointing in the direction of the MBH spin vector in the SSB frame
    S_ssb = [sin(ThetaKerr) * cos(PhiKerr), sin(ThetaKerr) * sin(PhiKerr), cos(ThetaKerr)]

    # compute unit vector pointing to observer (see text below Eq. 7 in https://arxiv.org/pdf/2104.04582)
    Theta_View = acos(dot(-R_ssb, S_ssb));
    Phi_View = -π/2;
    n_to_obs_src = [sin(Theta_View) * cos(Phi_View), sin(Theta_View) * sin(Phi_View), cos(Theta_View)];

    
    # calculate metric perturbations in source frame
    @inbounds Threads.@threads for t=1:nPoints
        @inbounds for i=1:3
            @inbounds for j=i:3
                hij[i, j][t] += 2.0 * Mij2[i, j][t] / r    # first term in Eq. 84 

                @inbounds for k=1:3
                    hij[i, j][t] += 2.0 * Mijk3[i, j, k][t] * n_to_obs_src[k] / (3.0r)    # second term in Eq. 84
    
                    @inbounds for l=1:3
                        hij[i, j][t] += 4.0 * (εkl[k, l][i] * Sij2[j, k][t] * n_to_obs_src[l] + εkl[k, l][j] * Sij2[i, k][t] * n_to_obs_src[l]) / (3.0r) + Mijkl4[i, j, k, l][t] * n_to_obs_src[k] * n_to_obs_src[l] / (6.0r)    # third and fourth terms in Eq. 84
        
                        @inbounds for m=1:3
                            hij[i, j][t] += (εkl[k, l][i] * Sijk3[j, k, m][t] * n_to_obs_src[l] * n_to_obs_src[m] + εkl[k, l][j] * Sijk3[i, k, m][t] * n_to_obs_src[l] * n_to_obs_src[m]) / (2.0r)
                        end
                    end
                end
            end
        end
    end

    hij[2, 1] = hij[1, 2]
    hij[3, 1] = hij[1, 3]
    hij[3, 2] = hij[2, 3]

    # compute tensor which projects into TT gauge (see text below Eq. 2 and Eq. 59 in https://arxiv.org/pdf/gr-qc/0202016)
    P = [δ(i, j) - n_to_obs_src[i] * n_to_obs_src[j] for i=1:3, j=1:3];
    Πijmn = [P[i, m] * P[j, n] - 0.5 * P[i,j] * P[m,n] for i=1:3, j=1:3, m=1:3, n=1:3];

    # compute wave polarization tensor (see Eqs. 22-23 in https://arxiv.org/pdf/1705.04259)
    p = [1.0, 0.0, 0.0] # x-axis
    q = [0, cos(Theta_View), sin(Theta_View)] # y-axis (after rotating source frame z-axis onto line of sight R)
    Hplus = [p[i] * p[j] - q[i] * q[j] for i=1:3, j=1:3]
    Hcross = [p[i] * q[j] + q[i] * p[j] for i=1:3, j=1:3]

    @inbounds for i = 1:3, j = 1:3
        for m=1:3, n=1:3
            hij_TT[i, j] += Πijmn[i, j, m, n] * hij[m, n]
        end

        hplus[:] += 0.5 * Hplus[i, j] * hij_TT[i, j]
        hcross[:] += 0.5 * Hcross[i, j] * hij_TT[i, j]
    end

    # normalize by mass ratio
    return hplus, hcross
end

function compute_wave_polarizations(nPoints::Int, r::Float64, ThetaSource::Float64, PhiSource::Float64, ThetaKerr::Float64, PhiKerr::Float64, Mij2::AbstractArray, mass_ratio::Float64)
    hij = [zeros(nPoints) for i=1:3, j=1:3];
    hij_TT = [zeros(nPoints) for i=1:3, j=1:3];
    hplus = zeros(nPoints);
    hcross = zeros(nPoints);
    
    # R_ssb ≡ unit vector pointing from solar system barycenter (SSB) in the direction of the EMRI system in the SSB frame
    R_ssb = [sin(ThetaSource) * cos(PhiSource), sin(ThetaSource) * sin(PhiSource), cos(ThetaSource)]
    
    # S_ssb ≡ unit vector pointing in the direction of the MBH spin vector in the SSB frame
    S_ssb = [sin(ThetaKerr) * cos(PhiKerr), sin(ThetaKerr) * sin(PhiKerr), cos(ThetaKerr)]

    # compute unit vector pointing to observer (see text below Eq. 7 in https://arxiv.org/pdf/2104.04582)
    Theta_View = acos(dot(-R_ssb, S_ssb));
    Phi_View = -π/2;
    n_to_obs_src = [sin(Theta_View) * cos(Phi_View), sin(Theta_View) * sin(Phi_View), cos(Theta_View)];

    
    # calculate metric perturbations in source frame
    @inbounds Threads.@threads for t=1:nPoints
        @inbounds for i=1:3
            @inbounds for j=i:3
                hij[i, j][t] = 2.0 * Mij2[i, j][t] / r    # first term in Eq. 84 
            end
        end
    end

    hij[2, 1] = hij[1, 2]
    hij[3, 1] = hij[1, 3]
    hij[3, 2] = hij[2, 3]

    # compute tensor which projects into TT gauge (see text below Eq. 2 and Eq. 59 in https://arxiv.org/pdf/gr-qc/0202016)
    P = [δ(i, j) - n_to_obs_src[i] * n_to_obs_src[j] for i=1:3, j=1:3];
    Πijmn = [P[i, m] * P[j, n] - 0.5 * P[i,j] * P[m,n] for i=1:3, j=1:3, m=1:3, n=1:3];

    # compute wave polarization tensor (see Eqs. 22-23 in https://arxiv.org/pdf/1705.04259)
    p = [1.0, 0.0, 0.0] # x-axis
    q = [0, cos(Theta_View), sin(Theta_View)] # y-axis (after rotating source frame z-axis onto line of sight R)
    Hplus = [p[i] * p[j] - q[i] * q[j] for i=1:3, j=1:3]
    Hcross = [p[i] * q[j] + q[i] * p[j] for i=1:3, j=1:3]

    @inbounds for i = 1:3, j = 1:3
        for m=1:3, n=1:3
            hij_TT[i, j] += Πijmn[i, j, m, n] * hij[m, n]
        end

        hplus[:] += 0.5 * Hplus[i, j] * hij_TT[i, j]
        hcross[:] += 0.5 * Hcross[i, j] * hij_TT[i, j]
    end

    # normalize by mass ratio
    return hplus, hcross
end

# returns plus and cross polarized waveforms with respect to the polarization vectors e_{ThetaObs}, e_{PhiObs} as defined in Eq. 21 in https://arxiv.org/pdf/gr-qc/0607007
@views function compute_wave_polarizations(nPoints::Int, r::Float64, ThetaObs::Float64, PhiObs::Float64, Mij2::AbstractArray, Mijk3::AbstractArray, Mijkl4::AbstractArray, Sij2::AbstractArray, Sijk3::AbstractArray, mass_ratio::Float64)
    hij = [zeros(nPoints) for i=1:3, j=1:3];
    hij_TT = [zeros(nPoints) for i=1:3, j=1:3];
    hplus = zeros(nPoints);
    hcross = zeros(nPoints);
    
    # n ≡ unit vector pointing in direction of far away observer
    nx = sin(ThetaObs) * cos(PhiObs)
    ny = sin(ThetaObs) * sin(PhiObs)
    nz = cos(ThetaObs)
    n_to_obs = [nx, ny, nz]
    
    # calculate perturbations (Eq. 84)
    @inbounds Threads.@threads for t=1:nPoints
        for i=1:3
            @inbounds for j=1:3
                hij[i, j][t] = 2.0 * Mij2[i, j][t] / r    # first term in Eq. 84 

                @inbounds for k=1:3
                    hij[i, j][t] += 2.0 * Mijk3[i, j, k][t] * n_to_obs[k] / (3.0r)    # second term in Eq. 84
    
                    @inbounds for l=1:3
                        hij[i, j][t] += 4.0 * (εkl[k, l][i] * Sij2[j, k][t] * n_to_obs[l] + εkl[k, l][j] * Sij2[i, k][t] * n_to_obs[l]) / (3.0r) + Mijkl4[i, j, k, l][t] * n_to_obs[k] * n_to_obs[l] / (6.0r)    # third and fourth terms in Eq. 84
        
                        @inbounds for m=1:3
                            hij[i, j][t] += (εkl[k, l][i] * Sijk3[j, k, m][t] * n_to_obs[l] * n_to_obs[m] + εkl[k, l][j] * Sijk3[i, k, m][t] * n_to_obs[l] * n_to_obs[m]) / (2.0r)
                        end
                    end
                end
            end
        end
        hplus[t] = Waveform.hplus(hij, ThetaObs, PhiObs, t)
        hcross[t] = Waveform.hcross(hij, ThetaObs, PhiObs, t)
    end
    
    # normalize by mass ratio
    return hplus, hcross
end


# returns plus and cross polarized waveforms with respect to the polarization vectors e_{ThetaObs}, e_{PhiObs} as defined in Eq. 21 in https://arxiv.org/pdf/gr-qc/0607007
@views function compute_wave_polarizations(nPoints::Int, r::Float64, ThetaObs::Float64, PhiObs::Float64, Mij2::AbstractArray, mass_ratio::Float64)
    hij = [zeros(nPoints) for i=1:3, j=1:3];
    hij_TT = [zeros(nPoints) for i=1:3, j=1:3];
    hplus = zeros(nPoints);
    hcross = zeros(nPoints);
    
    # n ≡ unit vector pointing in direction of far away observer
    nx = sin(ThetaObs) * cos(PhiObs)
    ny = sin(ThetaObs) * sin(PhiObs)
    nz = cos(ThetaObs)
    n_to_obs = [nx, ny, nz]

    # compute tensor which projects into TT gauge (see text below Eq. 2 and Eq. 59 in https://arxiv.org/pdf/gr-qc/0202016)
    P = [δ(i, j) - n_to_obs[i] * n_to_obs[j] for i=1:3, j=1:3];
    Πijmn = [P[i, m] * P[j, n] - 0.5 * P[i,j] * P[m,n] for i=1:3, j=1:3, m=1:3, n=1:3];

    # compute wave polarization tensor (see Eqs. 22-23 in https://arxiv.org/pdf/1705.04259)
    p = [1.0, 0.0, 0.0] # x-axis
    q = [0, cos(ThetaObs), sin(ThetaObs)] # y-axis (after rotating source frame z-axis onto line of sight R)
    Hplus = [p[i] * p[j] - q[i] * q[j] for i=1:3, j=1:3]
    Hcross = [p[i] * q[j] + q[i] * p[j] for i=1:3, j=1:3]
    
    # calculate perturbations (Eq. 84)
    @inbounds Threads.@threads for t=1:nPoints
        for i=1:3
            @inbounds for j=1:3
                hij[i, j][t] = 2.0 * Mij2[i, j][t] / r    # first term in Eq. 84 
            end
        end
        hplus[t] = Waveform.hplus(hij, ThetaObs, PhiObs, t)
        hcross[t] = Waveform.hcross(hij, ThetaObs, PhiObs, t)
    end
    
    # normalize by mass ratio
    return hplus, hcross
end

# project h into TT gauge (Reference: https://arxiv.org/pdf/gr-qc/0607007)
hΘΘ(h::AbstractArray, Θ::Float64, Φ::Float64, t::Int64) = (cos(Θ)^2) * (h[1, 1][t] * cos(Φ)^2 + h[1, 2][t] * sin(2Φ) + h[2, 2][t] * sin(Φ)^2) + h[3, 3][t] * sin(Θ)^2 - sin(2Θ) * (h[1, 3][t] * cos(Φ) + h[2, 3][t] * sin(Φ))   # Eq. 6.15

hΘΦ(h::AbstractArray, Θ::Float64, Φ::Float64, t::Int64) = cos(Θ) * (-0.5 * h[1, 1][t] * sin(2Φ) + h[1, 2][t] * cos(2Φ) + 0.5 * h[2, 2][t] * sin(2Φ)) + sin(Θ) * (h[1, 3][t] * sin(Φ) - h[2, 3][t] * cos(Φ))   # Eq. 6.16
hΦΦ(h::AbstractArray, Θ::Float64, Φ::Float64, t::Int64) = h[1, 1][t] * sin(Φ)^2 - h[1, 2][t] * sin(2Φ) + h[2, 2][t] * cos(Φ)^2   # Eq. 6.17

# define h_{+} and h_{×} components of GW (https://arxiv.org/pdf/gr-qc/0607007)
hplus(h::AbstractArray, Θ::Float64, Φ::Float64, t::Int64) = (1/2) *  (hΘΘ(h, Θ, Φ, t) - hΦΦ(h, Θ, Φ, t))
hcross(h::AbstractArray, Θ::Float64, Φ::Float64, t::Int64) = hΘΦ(h, Θ, Φ, t)


function h_plus_cross(hij::AbstractArray, Θ::Float64, Φ::Float64)
    nPoints = length(hij[1, 1])
    hplus = zeros(nPoints)
    hcross = zeros(nPoints)
    @inbounds Threads.@threads for i in 1:nPoints
        hplus[i] = Waveform.hplus(hij, Θ, Φ, i)
        hcross[i] = Waveform.hcross(hij, Θ, Φ, i)
    end
    return hplus, hcross
end

# Eq. 8 in https://arxiv.org/pdf/2104.04582
function rotate_to_SSB_frame(h_plus::Vector{Float64}, h_cross::Vector{Float64}, ThetaSource::Float64, PhiSource::Float64, ThetaK::Float64, PhiK::Float64)
    tan_psi_denominator = sin(ThetaK) * sin(PhiSource - PhiK)
    
    if abs(tan_psi_denominator) < 1e-10
        psi = π / 2
    else
        tan_psi_numerator = cos(ThetaSource) * sin(ThetaK) * cos(PhiSource - PhiK) - sin(ThetaSource) * cos(ThetaK)
        psi = -atan(tan_psi_numerator, tan_psi_denominator)
    end

    h_plus_SSB = h_plus * cos(2psi) - h_cross * sin(2psi)
    h_cross_SSB = h_plus * sin(2psi) + h_cross * cos(2psi)
    return h_plus_SSB, h_cross_SSB
end

function rotate_to_source_frame(h_plus::Vector{Float64}, h_cross::Vector{Float64}, ThetaSource::Float64, PhiSource::Float64, ThetaK::Float64, PhiK::Float64)
    tan_psi_denominator = sin(ThetaK) * sin(PhiSource - PhiK)
    
    if abs(tan_psi_denominator) < 1e-10
        psi = π / 2
    else
        tan_psi_numerator = cos(ThetaSource) * sin(ThetaK) * cos(PhiSource - PhiK) - sin(ThetaSource) * cos(ThetaK)
        psi = -atan(tan_psi_numerator, tan_psi_denominator)
    end

    h_plus_SSB = h_plus * cos(-2psi) - h_cross * sin(-2psi)
    h_cross_SSB = h_plus * sin(-2psi) + h_cross * cos(-2psi)
    return h_plus_SSB, h_cross_SSB
end

end