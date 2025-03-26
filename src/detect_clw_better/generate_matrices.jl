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

function generate_upwelling_readings(
    atm::AtmosphereParams,
    bands,
    clw_dist,
    t_sfc::Unitful.Temperature,
    emiss;
    res = Int(1e3),
)
    sfc = SurfaceParams(t_sfc, emiss, 1-emiss)
    down_rad = 0.0u"W/m^3"
    up_rad = 0.0u"W/m^3"
    
    cloud_atm = deepcopy(atm) 
    cloud_atm.clw .= clw_dist
    cloud_atm = upscale_atm_res(cloud_atm, res)

    function get_upwelling_rad(band, atm)
        model = RTEModel(
            sfc,
            atm,
            RadiativeParams(band, down_rad, up_rad, atm, 0.0u"°")
        )
        return solve_rte_only_TOA_SFC(model).up
    end

    cloudy_readings = get_upwelling_rad(bands, cloud_atm)
    return cloudy_readings 
end