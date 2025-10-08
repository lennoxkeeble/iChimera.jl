#=

    In this module we write functions which convert between the Chimera anomalies (ψ, χ, ϕ) to the AAK/AK anomalies (Φ, γ, alpha) as per Sec. IIIC https://arxiv.org/pdf/1705.04259

=#
module AnomalyConversion
using LinearAlgebra

function get_chi_phi(psi::Float64, gamma::Float64, alpha::Float64, θmin::Float64, sign_Lz::Int64, ThetaS::Float64, PhiS::Float64, ThetaK::Float64, PhiK::Float64)
    # compute S vector in solar system frame
    S = [sin(ThetaK) * cos(PhiK), sin(ThetaK) * sin(PhiK), cos(ThetaK)]
    R = [sin(ThetaS) * cos(PhiS), sin(ThetaS) * sin(PhiS), cos(ThetaS)]

    # compute L vector in solar system frame according to Eq. 12
    zhat = [0.0, 0.0, 1.0]
    inc = π/2 - sign_Lz * θmin
    inc < 1e-3 ? inc = 1e-3  : nothing
    L = S * cos(inc) + sin(inc) * (cos(alpha) * (zhat - dot(zhat, S) * S) / norm(zhat - dot(zhat, S) * S) + sin(alpha) * cross(S, zhat) / norm(cross(S, zhat)))

    # compute L-based coordinate frame basis vectors (where the components are with respect to the solar system frame) according to Eq. 47
    xhat_L = cross(L, S); xhat_L /= norm(xhat_L)
    yhat_L = dot(S, L) * L - S; yhat_L /= norm(yhat_L)
    zhat_L = L

    # compute S-based coordinate frame basis vectors (where the components are with respect to the solar system frame) according to Eq. 18
    xhat_S = cross(R, S); xhat_S /= norm(xhat_S)
    yhat_S = R - dot(R, S) * S; yhat_S /= norm(yhat_S)
    zhat_S = S

    # evaluate Eq. 48
    rL = [cos(psi + gamma), sin(psi + gamma), 0.0]

    QL = zeros(3, 3); QL[:, 1] = xhat_L; QL[:, 2] = yhat_L; QL[:, 3] = zhat_L
    QS = zeros(3, 3); QS[:, 1] = xhat_S; QS[:, 2] = yhat_S; QS[:, 3] = zhat_S

    rS = transpose(QS) * QL * rL # ≡ [sin(θ) cos(ϕ), sin(θ) sin(ϕ), cos(θ)]

    # compute BL coordinate theta at time zero
    cosθ0 = rS[3]; θ0 = acos(cosθ0)

    # compute chi0 from Eq. 35
    chi0 = acos(cosθ0 / cos(inc)) # here have set θmin -> inc since AAK built from Keplerian ellipses which will have constant inclination equal to θmin

    # compute phi0
    phi0 = asin(rS[2] / sin(θ0))
    
    return chi0, phi0
end

end