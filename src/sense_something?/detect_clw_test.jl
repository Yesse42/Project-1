cd(@__DIR__)
include("detect_clw.jl")

using Random

#Random.seed!(123)

band = 10u"GHz"

n_profiles = 500
errors = Float64[]

for _ in 1:n_profiles
    amount = rand((0.01:0.01:2)u"g/m^3")
    center = rand((1:0.1:3)u"km")
    spread = rand((0.1:0.1:1)u"km")
    rand_func(z) = amount * exp(-((z - center)/spread)^2)
    upwelling_radiation = get_upwelling_radiation(rand_func, band)
    clw_sensed = get_clw_from_10GHz(upwelling_radiation)
    clw_actual = trapz(standard_atm.z, rand_func.(standard_atm.z))
    push!(errors, abs((clw_sensed - clw_actual)/clw_actual))
end

mean_error = mean(errors)
println("Mean error in upwelling radiation: ", mean_error)
