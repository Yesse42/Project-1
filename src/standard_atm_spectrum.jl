cd(@__DIR__)
include("./test_cases/model_atm.jl")
include("external scripts/absorb_py.jl")

using Plots

visdir = "/Users/C837213770/Desktop/CSU HW/ATS 622/Project 1/vis"

standard_atm_first_level = first(eachrow(standard_atm))
T = standard_atm_first_level.T
P = standard_atm_first_level.Pressure
H20 = standard_atm_first_level.H20
clw =0.0u"g/m^3"

freq_range = range(1u"GHz", 300u"GHz"; length = 100_000)

k_abss = absorb_all.(freq_range, T, P, H20, clw)

plot(freq_range, k_abss, xlabel="Frequency", ylabel="Absorption Coefficient", title="Absorption Spectrum of Std. Atm. Layer 1", label = "", dpi = 500)
savefig(joinpath(visdir, "absorption_spectrum.png"))

