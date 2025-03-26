cd(@__DIR__)
include("../model_structures.jl")
include("../model_functions.jl")
include("../test_cases/model_atm.jl")

using Plots

first_standard_atm_layer = first(eachrow(standard_atm))
T = first_standard_atm_layer.T
P = first_standard_atm_layer.Pressure
#Get the frequency range
freqs = (1:0.25:300) * u"GHz"
#Get the cloud liquid water range
clws = (0:0.1:10) * u"g/m^3"
#Get the water vapor range
h2os = (0:0.1:10) * u"hPa"

#For the first layer of the standard atmosphere plot a heatmap showing the absorptivity as a function of frequency and cloud liquid water
clw_ks = absorb_all.(freqs', T, P, 0.0u"hPa", clws)
vapor_ks = absorb_all.(freqs', T, P, h2os, 0.0u"g/m^3")

# Display the heatmap for cloud liquid water absorptivity
heatmap(freqs, clws, clw_ks,
    xlabel="Frequency (GHz)", ylabel="Cloud Liquid Water (g/m^3)",
    title="Absorptivity vs Frequency and Cloud Liquid Water",
    colorbar_title="Absorptivity")
display(current())

# Display the heatmap for water vapor absorptivity
heatmap(freqs, h2os, vapor_ks,
    xlabel="Frequency (GHz)", ylabel="Water Vapor (hPa)",
    title="Absorptivity vs Frequency and Water Vapor",
    colorbar_title="Absorptivity")
display(current())