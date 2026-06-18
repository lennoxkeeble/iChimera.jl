using iChimera
using Test
using StaticArrays

function rr_derivative_fixture(; dxm_scale::Float64=0.0, rr_derivative_model=:partial_field)
    x = @MVector [1.2, -0.7, 0.45]
    dx = @MVector [0.08, -0.03, 0.05]
    d2x = @MVector [-0.004, 0.002, -0.001]

    Mij5 = @MMatrix [0.01 * (2i - j) for i = 1:3, j = 1:3]
    Mij6 = @MMatrix [0.02 * (i + 3j) for i = 1:3, j = 1:3]
    Mij7 = @MMatrix [-0.015 * (2i + j) for i = 1:3, j = 1:3]
    Mij8 = @MMatrix [0.005 * (i - 2j) for i = 1:3, j = 1:3]

    dxmMij5 = @MArray [dxm_scale * (i + 2j - k) for i = 1:3, j = 1:3, k = 1:3]
    dxmMij6 = @MArray [dxm_scale * (2i - j + 3k) for i = 1:3, j = 1:3, k = 1:3]
    dxmMij7 = @MArray [dxm_scale * (-i + 2j + k) for i = 1:3, j = 1:3, k = 1:3]

    Mijk7 = @MArray [0.003 * (i + 2j - k) for i = 1:3, j = 1:3, k = 1:3]
    Mijk8 = @MArray [-0.002 * (2i - j + k) for i = 1:3, j = 1:3, k = 1:3]
    dxmMijk7 = @MArray [dxm_scale * (i - 2j + 3k - l) for i = 1:3, j = 1:3, k = 1:3, l = 1:3]

    Sij5 = @MMatrix [0.004 * (i - j) for i = 1:3, j = 1:3]
    Sij6 = @MMatrix [-0.006 * (i + j) for i = 1:3, j = 1:3]
    dxmSij5 = @MArray [dxm_scale * (3i - j + k) for i = 1:3, j = 1:3, k = 1:3]

    Virr = @MVector zeros(3)
    dVrr_da = @MVector zeros(3)
    dVirr_dt = @MVector zeros(3)
    dVirr_da = @MMatrix zeros(3, 3)

    Vrr, dVrr_dt = iChimera.RRPotentials.compute_RR_potentials!(
        Virr,
        dVrr_da,
        dVirr_dt,
        dVirr_da,
        x,
        dx,
        d2x,
        Mij5,
        Mij6,
        Mij7,
        Mij8,
        dxmMij5,
        dxmMij6,
        dxmMij7,
        Mijk7,
        Mijk8,
        dxmMijk7,
        Sij5,
        Sij6,
        dxmSij5,
        0;
        rr_derivative_model=rr_derivative_model,
    )

    return (; Vrr, dVrr_dt, Virr=copy(Virr), dVrr_da=copy(dVrr_da), dVirr_dt=copy(dVirr_dt), dVirr_da=copy(dVirr_da))
end

@testset "RR derivative model normalization" begin
    @test iChimera.RRPotentials.normalize_rr_derivative_model(:partial_field) === :partial_field
    @test iChimera.RRPotentials.normalize_rr_derivative_model("corrected") === :partial_field
    @test iChimera.RRPotentials.normalize_rr_derivative_model(:legacy_worldline) === :legacy_worldline
    @test iChimera.RRPotentials.normalize_rr_derivative_model("legacy") === :legacy_worldline
    @test_throws ArgumentError iChimera.RRPotentials.normalize_rr_derivative_model(:not_a_derivative_model)
end

@testset "Partial-field RR derivatives ignore source-coordinate derivatives" begin
    partial_base = rr_derivative_fixture(dxm_scale=0.0)
    partial_shifted = rr_derivative_fixture(dxm_scale=0.75; rr_derivative_model=:partial_field)
    partial_alias = rr_derivative_fixture(dxm_scale=0.75; rr_derivative_model="corrected")

    @test partial_shifted.Vrr ≈ partial_base.Vrr
    @test partial_shifted.dVrr_dt ≈ partial_base.dVrr_dt
    @test partial_shifted.Virr ≈ partial_base.Virr
    @test partial_shifted.dVrr_da ≈ partial_base.dVrr_da
    @test partial_shifted.dVirr_dt ≈ partial_base.dVirr_dt
    @test partial_shifted.dVirr_da ≈ partial_base.dVirr_da
    @test partial_alias.dVrr_da ≈ partial_shifted.dVrr_da

    legacy_base = rr_derivative_fixture(dxm_scale=0.0; rr_derivative_model=:legacy_worldline)
    legacy_shifted = rr_derivative_fixture(dxm_scale=0.75; rr_derivative_model=:legacy_worldline)

    @test legacy_shifted.Vrr ≈ legacy_base.Vrr
    @test legacy_shifted.Virr ≈ legacy_base.Virr
    @test !isapprox(legacy_shifted.dVrr_da, legacy_base.dVrr_da)
    @test !isapprox(legacy_shifted.dVirr_da, legacy_base.dVirr_da)
end

@testset "Partial-field and legacy pulled-back derivatives are distinct" begin
    partial = rr_derivative_fixture(dxm_scale=0.4; rr_derivative_model=:partial_field)
    legacy = rr_derivative_fixture(dxm_scale=0.4; rr_derivative_model=:legacy_worldline)

    @test partial.Vrr ≈ legacy.Vrr
    @test partial.Virr ≈ legacy.Virr
    @test !isapprox(partial.dVrr_dt, legacy.dVrr_dt)
    @test !isapprox(partial.dVrr_da, legacy.dVrr_da)
    @test !isapprox(partial.dVirr_dt, legacy.dVirr_dt)
    @test !isapprox(partial.dVirr_da, legacy.dVirr_da)
end
