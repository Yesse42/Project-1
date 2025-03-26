include("model_structures.jl")
external_scripts_dir = "/Users/C837213770/Desktop/CSU HW/ATS 622/Project 1/src/external scripts"

using Unitful, DifferentialEquations, Interpolations, Plots, LinearAlgebra
import Unitful as U
import Base.Iterators as Itr

B(λ, T) = upreferred(2*U.h*U.c^2 / (λ^5 * (exp(U.h*U.c / (λ*U.k*T)) - 1)))

T_b(I, λ) = upreferred(U.h*U.c/(U.k*λ) * log((1+2*U.h*U.c^2/(I*λ^5)))^-1)

I_f(I0, B0, k0, Δs) = B0 - (B0 - I0)*exp(-k0*Δs)
new_rad_from_layer(B0, dτ) = B0 * (1 - exp(-dτ))

function solve_constant_Bs_ks(Δss, B_means, k_means, I0)
    Is = Vector{typeof(I0)}(undef, length(Δss) + 1)
    Is[1] = I0
    I = I0
    for (i, (Δs, B_mean, k_mean)) in enumerate(zip(Δss, B_means, k_means))
        I = I_f(I, B_mean, k_mean, Δs)
        Is[i + 1] = I
    end
    return Is
end

"Get the amount of radiation coming from each layer; as the idx of tau increases the light travels further along the path; note that tau is the accumulated optical depth along the path"
function get_radiation_from_layer(dτs, τs, B_means)
    I_from_layers = [new_rad_from_layer(B0, dτ) * exp(-(τs[end] - τs[i])) for (i, (dτ, B0)) in enumerate(zip(dτs, B_means))]
    return I_from_layers
end

function get_upwelling_weights(model::RTEModel, only_down = false)
    heights = model.atmosphere.z
    Δss = diff(heights)./cosd(model.radiative.view_zenith_angle)
    T_means = (model.atmosphere.T[1:end-1] .+ model.atmosphere.T[2:end]) ./ 2
    B_means = B.(model.radiative.λ, T_means)
    k_means = (model.radiative.k_abs[1:end-1] .+ model.radiative.k_abs[2:end]) ./ 2

    dτs = Δss .* k_means
    τs_up = cumsum(dτs)
    τs_down = cumsum(Itr.reverse(dτs))

    #First get the downwelling radiation
    down_radiation_amounts = reverse!(get_radiation_from_layer(Itr.reverse(dτs), τs_down, Itr.reverse(B_means)))

    final_I0_down = model.radiative.I0_down * exp(-τs_down[end])

    total_downwelling_rad = sum(down_radiation_amounts) + final_I0_down

    if only_down
        total_rad = sum(down_radiation_amounts) + final_I0_down
        return (layer_weights = upreferred.(down_radiation_amounts ./ total_rad),
                I0_down_weight = upreferred(final_I0_down / total_rad),
                total_rad = total_downwelling_rad)
    end

    #Now reflect the radiation from the surface, and attenuate through the whole atmosphere
    attenuated_from_sfc = exp(-τs_up[end])
    reflected_from_sfc_frac = model.surface.α_s * attenuated_from_sfc
    final_I0_down *= reflected_from_sfc_frac
    down_radiation_amounts .*= reflected_from_sfc_frac

    #Now get the upwelling radiation
    surface_emission = model.surface.ϵ_s * B(model.radiative.λ, model.surface.T_s) * attenuated_from_sfc
    final_I0_up = model.radiative.I0_up * attenuated_from_sfc
    up_radiation_amounts = get_radiation_from_layer(dτs, τs_up, B_means)
    rad_from_each_layer = up_radiation_amounts .+ down_radiation_amounts
    total_rad = sum(rad_from_each_layer) + final_I0_up + surface_emission + final_I0_down
    return (layer_weights = upreferred.(rad_from_each_layer ./ total_rad),
            sfc_emiss_weight = upreferred(surface_emission / total_rad),
            I0_up_weight = upreferred(final_I0_up / total_rad),
            I0_down_weight = upreferred(final_I0_down / total_rad),
            total_upwelling_rad = total_rad,
            total_downwelling_rad = total_downwelling_rad)
end

function solve_RTE(model::RTEModel)
    heights = model.atmosphere.z
    Δss = diff(heights)./cosd(model.radiative.view_zenith_angle)
    T_means = (model.atmosphere.T[1:end-1] .+ model.atmosphere.T[2:end]) ./ 2
    B_means = B.(model.radiative.λ, T_means)
    k_means = (model.radiative.k_abs[1:end-1] .+ model.radiative.k_abs[2:end]) ./ 2

    #Itr.reverse accounts for the downward direction
    down_Is = reverse!(solve_constant_Bs_ks(Itr.reverse(Δss), Itr.reverse(B_means), Itr.reverse(k_means), model.radiative.I0_down))

    #Now calulate the upward radiation
    I0_up = model.radiative.I0_up
    I0_up += model.surface.α_s * down_Is[begin]
    I0_up += model.surface.ϵ_s * B(model.radiative.λ, model.surface.T_s)
    up_Is = solve_constant_Bs_ks(Δss, B_means, k_means, I0_up)

    return (;down = down_Is, up = up_Is)
end

function solve_constant_Bs_ks_no_array(Δss, B_means, k_means, I0)
    I = I0
    for (Δs, B_mean, k_mean) in zip(Δss, B_means, k_means)
        I = I_f(I, B_mean, k_mean, Δs)
    end
    return I
end

function solve_rte_only_TOA_SFC(model::RTEModel)
    heights = model.atmosphere.z
    Δss = diff(heights)./cosd(model.radiative.view_zenith_angle)
    T_means = (model.atmosphere.T[1:end-1] .+ model.atmosphere.T[2:end]) ./ 2
    B_means = B.(model.radiative.λ, T_means)
    k_means = (model.radiative.k_abs[1:end-1] .+ model.radiative.k_abs[2:end]) ./ 2

    down_I = solve_constant_Bs_ks_no_array(Itr.reverse(Δss), Itr.reverse(B_means), Itr.reverse(k_means), model.radiative.I0_down)

    I0_up = model.radiative.I0_up
    I0_up += model.surface.α_s * down_I
    I0_up += model.surface.ϵ_s * B(model.radiative.λ, model.surface.T_s)

    up_I = solve_constant_Bs_ks_no_array(Δss, B_means, k_means, I0_up)
    return (;down = down_I, up = up_I)
end

function plot_RTE_sol(model, result; brightness = true, remove_units = true)
    heights = model.atmosphere.z
    λ = model.radiative.λ

    down_vals = up_vals = nothing
    xlabel = ylabel = ""

    if brightness
        down_vals = T_b.(result.down, λ)
        up_vals = T_b.(result.up, λ)
        ylabel = "Brightness Temperature (K)"
    else
        down_vals = uconvert.(U.W / U.m^2 / U.μm, result.down)
        up_vals = uconvert.(U.W / U.m^2 / U.μm, result.up)
        ylabel = "Radiance (W/m^2/μm)"
    end

    if remove_units
        heights = ustrip.(heights)
        down_vals = ustrip.(upreferred.(down_vals))
        up_vals = ustrip.(upreferred.(up_vals))
        xlabel = "Height (km)"
    else
        xlabel = "Height (km)"
    end

    p = plot(heights, down_vals, label = "Downwelling", xlabel = xlabel, ylabel = ylabel)
    plot!(p, heights, up_vals, label = "Upwelling")
    
    return p
end