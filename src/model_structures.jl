external_scripts_dir = "/Users/C837213770/Desktop/CSU HW/ATS 622/Project 1/src/external scripts"
include(joinpath(external_scripts_dir, "absorb_py.jl"))

using Interpolations

mutable struct SurfaceParams{T, S}
    T_s::T
    ϵ_s::S
    α_s::S
end

mutable struct AtmosphereParams{S2, S3, S4, S5, S6}
    z::Vector{S2}
    T::Vector{S3}
    P::Vector{S4}
    ρ_v::Vector{S5}
    clw::Vector{S6}
    function AtmosphereParams(z::Vector{S2}, T::Vector{S3}, P::Vector{S4}, ρ_v::Vector{S5}, clw::Vector{S6}) where {S2, S3, S4, S5, S6}
        @assert length(z) == length(T) == length(P) == length(ρ_v) == length(clw)
        @assert issorted(z)
        new{S2, S3, S4, S5, S6}(z, T, P, ρ_v, clw)
    end
end

function upscale_atm_res(atmosphere::AtmosphereParams, n_levels::Int)
    higher_res_z = collect(range(extrema(atmosphere.z)..., length=n_levels))
    fields = (:T, :P, :ρ_v, :clw)
    higher_res_fields = ntuple(field -> begin
        interp = LinearInterpolation(atmosphere.z, getfield(atmosphere, fields[field]), extrapolation_bc=Line())
        interp(higher_res_z)
    end, length(fields))

    return AtmosphereParams(
        higher_res_z,
        higher_res_fields...
    )
end

mutable struct RadiativeParams{S, S2, T1, T2, L, V}
    λ::S
    ν::S2
    I0_down::T1
    I0_up::T2
    k_abs::Vector{L}
    view_zenith_angle::V #degrees
    function RadiativeParams(λ::Unitful.Length, I0_down::T1, I0_up::T2, atmosphere::AtmosphereParams, view_zenith_angle::V) where {T1, T2, V}
        ν = Unitful.c / λ
        k_abs = compute_k_abs(ν, atmosphere)
        new{eltype(λ), typeof(ν), T1, T2, eltype(k_abs), V}(λ, ν, I0_down, I0_up, k_abs, view_zenith_angle)
    end
end

function RadiativeParams(ν::Unitful.Frequency, I0_down::T1, I0_up::T2, atmosphere::AtmosphereParams, view_zenith_angle::V) where {T1, T2, V}
    λ = Unitful.c / ν
    k_abs = compute_k_abs(ν, atmosphere)
    return RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)
end

function compute_k_abs(λ::S, atmosphere::AtmosphereParams) where {S}
    k_abs = absorb_all.(λ, atmosphere.T, atmosphere.P, atmosphere.ρ_v, atmosphere.clw)
    return k_abs
end

struct RTEModel{T1, T2, T3}
    surface::T1
    atmosphere::T2
    radiative::T3
end
