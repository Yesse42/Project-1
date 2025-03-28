cd(@__DIR__)
include("../model_structures.jl")
include("../model_functions.jl")
include("model_atm.jl")

visdir = "/Users/C837213770/Desktop/CSU HW/ATS 622/Project 1/vis"

using Plots, Unitful, Interpolations

# Define representative constants
myfreq = 183.0u"GHz"  # Frequency
λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
#Convert wavelength to frequency in GHz
I0_down = 0.0u"W/m^2/μm"  # Initial downwelling intensity
I0_up = 0.0u"W/m^2/μm" # Initial upwelling intensity
view_zenith_angle = 53u"°"  # Zenith angle, degrees

clw = fill(0.0u"g/m^3", length(standard_atm.Hgt))

standard_atm = AtmosphereParams(
    standard_atm.Hgt,
    standard_atm.T,
    standard_atm.Pressure,
    standard_atm.H20,
    clw
)

# upscale the atmosphere resolution
atmosphere = upscale_atm_res(standard_atm, Int(1e5))

# Create the surface parameters
emiss = 0.5
surface = SurfaceParams(first(standard_atm.T), emiss, 1-emiss)

# Create the radiative parameters
radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)

# Create the RTE model
model = RTEModel(surface, atmosphere, radiative)

# Solve the RTE
solution = solve_RTE(model)

#convert solution to brightness temperature
solution = (down = T_b.(solution.down, model.radiative.λ),
             up = T_b.(solution.up, model.radiative.λ))
# Plot the results
domain = model.atmosphere.z
plot(domain, solution.down, label="Downwelling Intensity")
plot!(domain, solution.up, label="Upwelling Intensity"; ylims = (220, 290))
title!("$(myfreq) RTE, Standard Atm, θ = $(view_zenith_angle)")
xlabel!("Height (m)")
ylabel!("Brightness Temperature (K)")
display(plot!(; legend = :topright, dpi = 500))
savefig(joinpath(visdir, "RTE_solution_183_GHz.png"))

quicksol = solve_rte_only_TOA_SFC(model)
display(T_b(quicksol.up, λ))
display(T_b(quicksol.down, λ))
