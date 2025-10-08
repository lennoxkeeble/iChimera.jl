#=

    In this module we compute derivatives of Mino time, λ, with respect to coordinate time t, e.g., λ^{(n)}(t), for n = 1, 2, 3, 4, 5, 6, 7, 8, using t^{(n)}(λ) via the chain rule

=#
module MinoTimeDerivs
dλ_dt(dt_dλ::Float64) = 1/dt_dλ

d2λ_dt(dt_dλ::Float64, d2t_dλ::Float64) = -(d2t_dλ/dt_dλ^3)

d3λ_dt(dt_dλ::Float64, d2t_dλ::Float64, d3t_dλ::Float64) = (3*d2t_dλ^2 - dt_dλ*d3t_dλ)/dt_dλ^5

d4λ_dt(dt_dλ::Float64, d2t_dλ::Float64, d3t_dλ::Float64, d4t_dλ::Float64) = -((15*d2t_dλ^3 - 10*dt_dλ*d2t_dλ*d3t_dλ + dt_dλ^2*d4t_dλ)/dt_dλ^7)

d5λ_dt(dt_dλ::Float64, d2t_dλ::Float64, d3t_dλ::Float64, d4t_dλ::Float64, d5t_dλ::Float64) = (105*d2t_dλ^4 - 105*dt_dλ*d2t_dλ^2*d3t_dλ +
5*dt_dλ^2*(2*d3t_dλ^2 + 3*d2t_dλ*d4t_dλ) - dt_dλ^3*d5t_dλ)/dt_dλ^9

d6λ_dt(dt_dλ::Float64, d2t_dλ::Float64, d3t_dλ::Float64, d4t_dλ::Float64, d5t_dλ::Float64, d6t_dλ::Float64) = (7*(-135*d2t_dλ^5 + 180*dt_dλ*d2t_dλ^3*d3t_dλ -
30*dt_dλ^2*d2t_dλ^2*d4t_dλ + 5*dt_dλ^3*d3t_dλ*d4t_dλ + dt_dλ^2*d2t_dλ*(-40*d3t_dλ^2 + 3*dt_dλ*d5t_dλ)) - dt_dλ^4*d6t_dλ)/dt_dλ^11

d7λ_dt(dt_dλ::Float64, d2t_dλ::Float64, d3t_dλ::Float64, d4t_dλ::Float64, d5t_dλ::Float64, d6t_dλ::Float64, d7t_dλ::Float64) = (7*(1485*d2t_dλ^6 - 2475*dt_dλ*d2t_dλ^4*d3t_dλ + 
450*dt_dλ^2*d2t_dλ^3*d4t_dλ + 18*dt_dλ^2*d2t_dλ^2*(50*d3t_dλ^2 - 3*dt_dλ*d5t_dλ) + dt_dλ^3*(-40*d3t_dλ^3 + 5*dt_dλ*d4t_dλ^2 +
8*dt_dλ*d3t_dλ*d5t_dλ) + 4*dt_dλ^3*d2t_dλ*(-45*d3t_dλ*d4t_dλ + dt_dλ*d6t_dλ)) - dt_dλ^5*d7t_dλ)/dt_dλ^13

d8λ_dt(dt_dλ::Float64, d2t_dλ::Float64, d3t_dλ::Float64, d4t_dλ::Float64, d5t_dλ::Float64, d6t_dλ::Float64, d7t_dλ::Float64, d8t_dλ::Float64) = (-135135*d2t_dλ^7 + 270270*dt_dλ*d2t_dλ^5*d3t_dλ -
51975*dt_dλ^2*d2t_dλ^4*d4t_dλ + 6930*dt_dλ^2*d2t_dλ^3*(-20*d3t_dλ^2 + dt_dλ*d5t_dλ) - 630*dt_dλ^3*d2t_dλ^2*(-55*d3t_dλ*d4t_dλ +
dt_dλ*d6t_dλ) + dt_dλ^3*d2t_dλ*(-35*(-440*d3t_dλ^3 + 45*dt_dλ*d4t_dλ^2 + 72*dt_dλ*d3t_dλ*d5t_dλ) + 36*dt_dλ^2*d7t_dλ) + dt_dλ^4*(-2100*d3t_dλ^2*d4t_dλ + 42*dt_dλ*(3*d4t_dλ*d5t_dλ +
2*d3t_dλ*d6t_dλ) - dt_dλ^2*d8t_dλ))/dt_dλ^15

d9λ_dt(dt_dλ::Float64, d2t_dλ::Float64, d3t_dλ::Float64, d4t_dλ::Float64, d5t_dλ::Float64, d6t_dλ::Float64, d7t_dλ::Float64, d8t_dλ::Float64, d9t_dλ::Float64) =(2027025*d2t_dλ^8 - 4729725*dt_dλ*d2t_dλ^6*d3t_dλ + 945945*dt_dλ^2*d2t_dλ^5*d4t_dλ - 45045*dt_dλ^2*d2t_dλ^4*(-70*d3t_dλ^2 + 3*dt_dλ*d5t_dλ) + 13860*dt_dλ^3*d2t_dλ^3*(-65*d3t_dλ*d4t_dλ + dt_dλ*d6t_dλ) - 165*dt_dλ^3*d2t_dλ^2*(3640*d3t_dλ^3 - 315*dt_dλ*d4t_dλ^2 - 504*dt_dλ*d3t_dλ*d5t_dλ + 6*dt_dλ^2*d7t_dλ) + 15*dt_dλ^4*d2t_dλ*(462*d4t_dλ*(20*d3t_dλ^2 - dt_dλ*d5t_dλ) - 308*dt_dλ*d3t_dλ*d6t_dλ + 3*dt_dλ^2*d8t_dλ) + dt_dλ^4*(15400*d3t_dλ^4 - 4620*dt_dλ*d3t_dλ^2*d5t_dλ + 15*dt_dλ*d3t_dλ*(-385*d4t_dλ^2 + 8*dt_dλ*d7t_dλ) + dt_dλ^2*(126*d5t_dλ^2 + 210*d4t_dλ*d6t_dλ - dt_dλ*d9t_dλ)))/dt_dλ^17

end

MMA_vals = [-2.6426974176218043, -2.5784743561419488, -31.674992783332012,
-307.2930469891393, -5287.8001371780865, -106520.6074845081,
-2.707061276046853*1e6, -8.029580083179496*1e7,
-2.758801010302011*1e9]

Dλ(n) = sin(n+1) * cos(n+1)
Df(n) = sin(5n+2) * cos(5n+2)

dt_dλ = Dλ(1)
d2t_dλ = Dλ(2)
d3t_dλ = Dλ(3)
d4t_dλ = Dλ(4)
d5t_dλ = Dλ(5)
d6t_dλ = Dλ(6)
d7t_dλ = Dλ(7)
d8t_dλ = Dλ(8)
d9t_dλ = Dλ(9)

MinoTimeDerivs.dλ_dt(dt_dλ) ≈ MMA_vals[1]
MinoTimeDerivs.d2λ_dt(dt_dλ, d2t_dλ) ≈ MMA_vals[2]
MinoTimeDerivs.d3λ_dt(dt_dλ, d2t_dλ, d3t_dλ) ≈ MMA_vals[3]
MinoTimeDerivs.d4λ_dt(dt_dλ, d2t_dλ, d3t_dλ, d4t_dλ) ≈ MMA_vals[4]
MinoTimeDerivs.d5λ_dt(dt_dλ, d2t_dλ, d3t_dλ, d4t_dλ, d5t_dλ) ≈ MMA_vals[5]
MinoTimeDerivs.d6λ_dt(dt_dλ, d2t_dλ, d3t_dλ, d4t_dλ, d5t_dλ, d6t_dλ) ≈ MMA_vals[6]
MinoTimeDerivs.d7λ_dt(dt_dλ, d2t_dλ, d3t_dλ, d4t_dλ, d5t_dλ, d6t_dλ, d7t_dλ) ≈ MMA_vals[7]
MinoTimeDerivs.d8λ_dt(dt_dλ, d2t_dλ, d3t_dλ, d4t_dλ, d5t_dλ, d6t_dλ, d7t_dλ, d8t_dλ) ≈ MMA_vals[8]
MinoTimeDerivs.d9λ_dt(dt_dλ, d2t_dλ, d3t_dλ, d4t_dλ, d5t_dλ, d6t_dλ, d7t_dλ, d8t_dλ, d9t_dλ) ≈ MMA_vals[9]