cd(@__DIR__)
include("../model_structures.jl")
include("../model_functions.jl")
include("model_atm.jl")

testvisdir = "../../vis/tests"

using Polynomials, Plots

clw = fill(0.0u"g/m^3", length(standard_atm.Hgt))

#First get the standard atmosphere
standard_atm = AtmosphereParams(
    standard_atm.Hgt,
    standard_atm.T,
    standard_atm.Pressure,
    standard_atm.H20,
    clw
)

no_air_model = 
let
    atm = deepcopy(standard_atm)
    atm.ρ_v .= 0.0u"hPa"
    atm.P .= 0u"hPa"
    myfreq = 183.0u"GHz"  # Frequency
    λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
    #Convert wavelength to frequency in GHz
    I0_down = B(λ, 280u"K")  # Initial downwelling intensity
    I0_up = zero(typeof(I0_down)) # Initial upwelling intensity
    view_zenith_angle = 0u"°"  # Zenith angle, degrees

    atmosphere = upscale_atm_res(atm, Int(1e3))

    # Create the surface parameters
    emiss = 1.
    surface = SurfaceParams(first(standard_atm.T), emiss, 1-emiss)

    # Create the radiative parameters
    radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)

    # Create the RTE model
    model = RTEModel(surface, atmosphere, radiative)

    #Solve the RTE, and plot
    solution = solve_RTE(model)
    println("Surface Temperature, no atm: ", surface.T_s)
    plot = plot_RTE_sol(model, solution)
    plot!(title="No Air Model")
    savefig(joinpath(testvisdir, "no_air_model.png"))
end

water_vapor_discontinuity_model = 
let 
    atm = deepcopy(standard_atm)
    atm.ρ_v .= 0.0u"hPa"
    myfreq = 183.0u"GHz"  # Frequency
    λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
    #Convert wavelength to frequency in GHz
    I0_down = 0.0u"W/m^2/μm"  # Initial downwelling intensity
    I0_up = 0.0u"W/m^2/μm" # Initial upwelling intensity
    view_zenith_angle = 0u"°"  # Zenith angle, degrees

    atmosphere = upscale_atm_res(standard_atm, Int(1e3))
    atmosphere.ρ_v[isapprox.(atmosphere.z, 4u"km"; rtol = 0, atol = 1.0u"km")] .= 100u"hPa"

    # Create the surface parameters
    emiss = 0.5
    surface = SurfaceParams(first(standard_atm.T), emiss, 1-emiss)

    # Create the radiative parameters
    radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)

    # Create the RTE model
    model = RTEModel(surface, atmosphere, radiative)

    #Solve, plot, and then add the model temperature and H20 profiles to the plot
    solution = solve_RTE(model)
    plot = plot_RTE_sol(model, solution)
    plot!(plot, atmosphere.z .|> ustrip, atmosphere.T .|> ustrip, label="Temperature Profile")
    plot!(title="Water Vapor Discontinuity Model")
    savefig(joinpath(testvisdir, "water_vapor_discontinuity_model.png"))
end

full_wv_model = 
let 
    atm = deepcopy(standard_atm)
    atm.ρ_v .= 20.0u"hPa"
    myfreq = 183.0u"GHz"  # Frequency
    λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
    #Convert wavelength to frequency in GHz
    I0_down = B(λ, 280u"K")  # Initial downwelling intensity
    I0_up = zero(typeof(I0_down)) # Initial upwelling intensity
    view_zenith_angle = 0u"°"  # Zenith angle, degrees

    atmosphere = upscale_atm_res(atm, Int(1e3))

    # Create the surface parameters
    emiss = 1.
    surface = SurfaceParams(first(standard_atm.T), emiss, 1-emiss)

    # Create the radiative parameters
    radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)

    # Create the RTE model
    model = RTEModel(surface, atmosphere, radiative)

    #Solve the RTE, and save the plot
    solution = solve_RTE(model)
    plot = plot_RTE_sol(model, solution)
    plot!(title="Full Water Vapor Model")
    plot!(plot, atmosphere.z .|> ustrip, atmosphere.T .|> ustrip, label="Temperature Profile")
    savefig(joinpath(testvisdir, "full_wv_model.png"))
end

high_cloud_model = 
let 
    atm = deepcopy(standard_atm)
    atm.ρ_v .= 0.0u"hPa"
    myfreq = 183.0u"GHz"  # Frequency
    λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
    #Convert wavelength to frequency in GHz
    I0_down = 0.0u"W/m^2/μm"  # Initial downwelling intensity
    I0_up = 0.0u"W/m^2/μm" # Initial upwelling intensity
    view_zenith_angle = 0u"°"  # Zenith angle, degrees

    atmosphere = upscale_atm_res(atm, Int(1e3))
    atmosphere.clw[isapprox.(atmosphere.z, 5u"km"; rtol = 0, atol = 0.5u"km")] .= 0.1u"g/m^3"

    # Create the surface parameters
    emiss = 1.
    surface = SurfaceParams(first(standard_atm.T), emiss, 1-emiss)

    # Create the radiative parameters
    radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)

    # Create the RTE model
    model = RTEModel(surface, atmosphere, radiative)

    #Solve, plot, and then add the model temperature and cloud profiles to the plot
    solution = solve_RTE(model)
    plot = plot_RTE_sol(model, solution)
    plot!(title="High Cloud Model")
    savefig(joinpath(testvisdir, "high_cloud_model.png"))
end

no_wv_model = 
let
    atm = deepcopy(standard_atm)
    atm.ρ_v .= 0.0u"hPa"
    myfreq = 183.0u"GHz"  # Frequency
    λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
    #Convert wavelength to frequency in GHz
    I0_down = B(λ, 280u"K")  # Initial downwelling intensity
    I0_up = zero(typeof(I0_down)) # Initial upwelling intensity
    view_zenith_angle = 0u"°"  # Zenith angle, degrees

    atmosphere = upscale_atm_res(atm, Int(1e3))

    # Create the surface parameters
    emiss = 1.
    surface = SurfaceParams(first(standard_atm.T), emiss, 1-emiss)

    # Create the radiative parameters
    radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)

    # Create the RTE model
    model = RTEModel(surface, atmosphere, radiative)

    #Solve the RTE, and save the plot
    solution = solve_RTE(model)
    println("Surface Temperature, no water vapor: ", surface.T_s)
    plot = plot_RTE_sol(model, solution)
    plot!(title="No Water Vapor Model")
    savefig(joinpath(testvisdir, "no_wv_model.png"))
end

uniform_wv_model = 
let 
    atm = deepcopy(standard_atm)
    atm.T .= 280.0u"K"  # Constant temperature
    atm.P .= 1013.25u"hPa"  # Constant pressure
    myfreq = 183.0u"GHz"  # Frequency
    λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
    view_zenith_angle = 0u"°"  # Zenith angle, degrees
    atm.ρ_v .= 0.1u"hPa"

    atmosphere = upscale_atm_res(atm, Int(1e3))  # Uniform water vapor profile

    #Now make the temperature vary sinusoidally
    wavenum = 3
    T0 = 280.0u"K"
    dT = 10.0u"K"
    wavenum = 3 * 2pi / maximum(atmosphere.z)
    atmosphere.T .= T0 .+ dT .* sin.(wavenum * atmosphere.z)

    I0_down = B(λ, atmosphere.T[end])  # Initial downwelling intensity
    I0_up = zero(typeof(I0_down)) # Initial upwelling intensity

    # Create the surface parameters
    surface = SurfaceParams(first(standard_atm.T), 0.0, 1.0)

    # Create the radiative parameters
    radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)

    # Create the RTE model
    model = RTEModel(surface, atmosphere, radiative)

    #Solve the RTE, and save the plot
    solution = solve_RTE(model)
    plot = plot_RTE_sol(model, solution)
    plot!(title="Sinusoidal T Model")
    plot!(atmosphere.z .|> ustrip, atmosphere.T .|> ustrip, label="T Profile")
    savefig(joinpath(testvisdir, "sinusoidal_t_model.png"))

    println("sinusoid")
    println("k_abs: ", radiative.k_abs[1])
    println("Wavenumber of the sin function: ", wavenum * 2pi / maximum(atmosphere.z))
    println("Amplitude of the sin function: ", 10.0u"K")
end

let
    zenith_angles = 0:1:80

    myfreq = 180u"GHz"  # Frequency
    λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
    I0_down = 1.0u"W/m^2/μm"  # Initial downwelling intensity
    I0_up = zero(typeof(I0_down)) # Initial upwelling intensity

    downwelling_radiation = typeof(I0_down)[]

    atm = deepcopy(standard_atm)
    actual_optical_depth = nothing

    # Create the surface parameters
    emiss = 1
    surface = SurfaceParams(first(standard_atm.T), emiss, 1-emiss)

    for angle in zenith_angles

        view_zenith_angle = angle * u"°"  # Zenith angle, degrees

        atmosphere = upscale_atm_res(atm, Int(1e3))

        # Create the radiative parameters
        radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)
        actual_optical_depth = diff(atmosphere.z) ⋅ (radiative.k_abs[1:end-1] .+ radiative.k_abs[2:end]) ./ 2
        
        #Now set the atmosphere temp to 1K everywhere for no emission
        atmosphere.T .= 0.01u"K"

        # Create the RTE model
        model = RTEModel(surface, atmosphere, radiative)

        # Solve the RTE
        solution = solve_RTE(model)
        
        # Extract the downwelling radiation at the ground
        push!(downwelling_radiation, solution.down[1])
    end

    downwelling_radiation = uconvert.(u"W/m^2/μm", downwelling_radiation)
    intensity_units = unit(downwelling_radiation[1])
    seczenith = secd.(zenith_angles)

    linear_fit = fit(seczenith, log.(ustrip.(downwelling_radiation)), 1).coeffs
    τ = -linear_fit[2]
    println("τ = ", τ)
    println("Actual Optical Depth: ", actual_optical_depth)
    I0 = exp(linear_fit[1]) * intensity_units
    println("I0 = ", I0)
    println("Original Downwelling Radiation: ", uconvert(u"W/m^2/μm", I0_down))

    p = plot(seczenith, log.(ustrip.(downwelling_radiation)), label="Downwelling Radiation", xlabel="Sec Zenith Angle", ylabel="Log Downwelling Radiation (W/m^2/μm)", title="Downwelling Radiation vs Zenith Angle")
    plot!(seczenith, linear_fit[1] .+ linear_fit[2] .* seczenith, label="Linear Fit")
    savefig(joinpath(testvisdir, "downwelling_radiation_vs_zenith_angle.png"))
end

let
    frequencies = (1.0:0.01:300.0)u"GHz"  # Frequencies in GHz
    atm = upscale_atm_res(standard_atm, Int(200))
    function calculate_upwelling_brightness_temp(myfreq)  # Frequency
        λ = Unitful.c / myfreq  # Convert wavelength to frequency in GHz
        I0_down = B(λ, 0.00001u"K")  # Initial downwelling intensity
        I0_up = zero(typeof(I0_down))  # Initial upwelling intensity
        view_zenith_angle = 0u"°"  # Zenith angle, degrees

        atmosphere = atm

        # Create the surface parameters
        emiss = 1.0 # Surface emissivity
        surface = SurfaceParams(first(standard_atm.T), emiss, 1-emiss)

        # Create the radiative parameters
        radiative = RadiativeParams(λ, I0_down, I0_up, atmosphere, view_zenith_angle)

        # Create the RTE model
        model = RTEModel(surface, atmosphere, radiative)

        # Solve the RTE and calculate upwelling brightness temperature
        solution = solve_rte_only_TOA_SFC(model).up
        brightness_temp = T_b(solution, λ)
        return ustrip(brightness_temp)  # Return the result
    end

    upwelling_brightness_temps = calculate_upwelling_brightness_temp.(frequencies)

    # Plot the upwelling brightness temperature vs frequency
    p = plot(
        frequencies,
        collect(upwelling_brightness_temps),
        xlabel="Frequency",
        ylabel="Upwelling Brightness Temperature (K)",
        title="Upwelling Brightness Temperature vs Frequency",
        label="Brightness Temperature",
        dpi = 500
    )
    savefig(joinpath(testvisdir, "upwelling_brightness_temperature_vs_frequency.png"))
end