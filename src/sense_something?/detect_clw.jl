cd(@__DIR__)
include("../model_structures.jl")
include("../model_functions.jl")
include("../test_cases/model_atm.jl")

using Plots, DataFrames, BenchmarkTools, Trapz, Polynomials

visdir = joinpath(@__DIR__, "../..", "vis", "clw")

clw = fill(0.0u"g/m^3", length(standard_atm.Hgt))

standard_atm = AtmosphereParams(
    standard_atm.Hgt,
    standard_atm.T,
    standard_atm.Pressure,
    standard_atm.H20,
    clw
)

standard_atm = upscale_atm_res(standard_atm, 1000)

#Get some bands
bands = [89, 30, 10] * u"GHz"

"Make a function to get the upwelling radiation in the standard atmosphere for a given clw function"
function get_upwelling_radiation(clw_func, band::Unitful.Frequency; atm = standard_atm, t_sfc = 288u"K", emiss = 0.5)
    clw = clw_func.(atm.z)
    atm.clw .= clw
    rad_units = u"W/m^3"
    model = RTEModel(
        SurfaceParams(t_sfc, emiss, 1-emiss),
        standard_atm,
        RadiativeParams(band, 0.0*rad_units, 0.0*rad_units, standard_atm, 0.0u"°")
    )   
    return solve_rte_only_TOA_SFC(model).up
end

concs = (0.05:0.1:3) .* u"g/m^3"

cloudfuncs = reduce(vcat, [[(z -> if z<= 1000u"m" conc else 0.0u"g/m^3" end), (z-> if (3000u"m" <= z <= 4000u"m") conc else 0.0u"g/m^3" end), (z-> conc * exp(-((z - 2000u"m")/1u"km")^2))] for conc in concs])

colors = reduce(vcat, [[:red, :green, :blue] for _ in concs])
color_names = ["Low Cloud", "High Cloud", "Mid Level Gaussian Cloud"]

#Calculate the integrated water for each cloudfunc
integrated_water = [upreferred(trapz(standard_atm.z, func.(standard_atm.z))) for func in cloudfuncs]

#Get the upwelling radiation for each cloudfunc
upwelling_radiation = [collect(get_upwelling_radiation.(cloudfuncs, band)) for band in bands]

for (i, band) in enumerate(bands)
    p = scatter(
        uconvert.(u"W/m^2/μm", upwelling_radiation[i]),
        integrated_water,
        c = colors,
        label = "",
        xlabel = "Upwelling Radiation",
        ylabel = "Integrated Water",
        title = "Integrated Water vs Upwelling Radiation for $band",
        xlims = (-Inf, 1.6e-12),
        ylims = (-Inf, 2.4),
        dpi = 500
    )    #Add a legend
    for (j, color) in enumerate([:red, :green, :blue])
        scatter!(p, [NaN], [NaN], label = color_names[j], c = color)
    end
    display(p)
    if i == 3
        savefig(p, joinpath(visdir, "clw_vs_upwelling_10GHz.png"))
    end
end

best_idx = 3
lin_fit = coeffs(fit(Polynomial,
    ustrip.(upwelling_radiation[best_idx][:]),
    ustrip.(integrated_water[:]),
    1
))
const intercept = lin_fit[1] * unit(integrated_water[1])
const slope = lin_fit[2] * (unit(upwelling_radiation[best_idx][1]) / unit(integrated_water[1]))^-1
println("Intercept: $intercept")
println("Slope: $slope")

function get_clw_from_10GHz(reading)
    return intercept + slope * reading
end