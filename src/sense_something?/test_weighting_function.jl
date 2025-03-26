cd(@__DIR__)
include("../model_structures.jl")
include("../model_functions.jl")
include("../test_cases/model_atm.jl")

using Plots, DataFrames, BenchmarkTools

visdir = joinpath(@__DIR__, "../..", "vis", "weighting")

clw = fill(0.0u"g/m^3", length(standard_atm.Hgt))

# Function to create atmosphere with modified water vapor
function create_atmosphere(water_vapor_factor)
    AtmosphereParams(
        standard_atm.Hgt,
        standard_atm.T,
        standard_atm.Pressure,
        standard_atm.H20 .* water_vapor_factor,
        clw
    )
end

# Function to get weighting functions
function get_weights(band; atm)
    λ = Unitful.c / band  # Convert wavelength to frequency in GHz
    I0_down = 0.0u"W/m^3" # Initial downwelling intensity
    I0_up = I0_down # Initial upwelling intensity
    view_zenith_angle = 0u"°"  # Zenith angle, degrees

    emiss = 0.5
    surface = SurfaceParams(280.0u"K", emiss, 1-emiss)

    radiative = RadiativeParams(λ, I0_down, I0_up, atm, view_zenith_angle)

    model = RTEModel(surface, atm, radiative)

    sol = get_upwelling_weights(model)
    return sol
end

# Define parameters
nbands = 15
start_freq = 170u"GHz"
furthest_freq = 183.31u"GHz"
bands = collect(range(start_freq, stop=furthest_freq, length=nbands))

#push!(bands, 23.8u"GHz")
#push!(bands, 89u"GHz")

# Define water vapor factors and corresponding titles
water_vapor_factors = [1.0, 0.5, 2.0]
titles = [
    "Upwelling Radiation Weighting Functions (Standard Atmosphere)",
    "Upwelling Radiation Weighting Functions (Halved Std. Water Vapor)",
    "Upwelling Radiation Weighting Functions (Doubled Std. Water Vapor)"
]
filenames = [
    "weighting_functions_normal.png",
    "weighting_functions_halved.png",
    "weighting_functions_doubled.png"
]

n_layers = Int(1e4)

standard = upscale_atm_res(create_atmosphere(1), n_layers)

zvals = (standard.z[1:end-1] .+ standard.z[2:end]) ./ 2

lowest_20_km = findall(0u"m" .< zvals .<= 20_000u"m")

dz = standard.z[2] - standard.z[1]

# Loop through water vapor factors
for (factor, title, filename) in zip(water_vapor_factors, titles, filenames)
    # Create atmosphere
    atm = upscale_atm_res(create_atmosphere(factor), n_layers)
    
    # Compute weights
    weights = get_weights.(bands, atm=atm)
    
    # Plot
    plot(
        zvals[lowest_20_km], 
        reduce(hcat, [weight.layer_weights for weight in weights])[lowest_20_km, :] ./ dz, 
        xlabel="Height", 
        ylabel="Density", 
        label = "", 
        linez = permutedims(bands), 
        title=title, 
        cmap = :viridis, 
        colorbar_title="Frequency (GHz)",
        titlefontsize = 10,
        dpi = 500
    )
    #Also print out a dataframe of the sfc emission fraction for each band
    df = DataFrame(
        Frequency = bands,
        Surface_Emission_Fraction = [weight.sfc_emiss_weight for weight in weights]
    )
    println(factor)
    display(df)

    savefig(joinpath(visdir, filename))
end

# Define a function to remove water vapor above 1km
function remove_water_vapor_above_1km(atm)
    modified_h2o = atm.ρ_v
    modified_h2o[atm.z .> 1000u"m"] .= 0.0u"hPa"
    AtmosphereParams(atm.z, atm.T, atm.P, modified_h2o, atm.clw)
end

# Create atmosphere with no water vapor above 1km
atm_no_h2o_above_1km = remove_water_vapor_above_1km(standard)

# Compute weights for the modified atmosphere
weights_no_h2o_above_1km = get_weights.(bands, atm=atm_no_h2o_above_1km)

lowest_2_km = findall(0u"m" .< zvals .<= 2_000u"m")

# Plot for the modified atmosphere
plot(
    zvals[lowest_2_km], 
    reduce(hcat, [weight.layer_weights for weight in weights_no_h2o_above_1km])[lowest_2_km, :] ./ dz, 
    xlabel="Height", 
    ylabel="Density", 
    label = "", 
    linez = permutedims(bands), 
    title="Upwelling Radiation Weighting Functions (No Water Vapor Above 1km)", 
    cmap = :viridis, 
    colorbar_title="Frequency (GHz)",
    titlefontsize = 10,
    dpi = 500
)

# Print out a dataframe of the surface emission fraction for each band
df_no_h2o_above_1km = DataFrame(
    Frequency = bands,
    Surface_Emission_Fraction = [weight.sfc_emiss_weight for weight in weights_no_h2o_above_1km]
)
println("No Water Vapor Above 1km")
display(df_no_h2o_above_1km)

savefig(joinpath(visdir, "weighting_functions_no_h2o_above_1km.png"))


