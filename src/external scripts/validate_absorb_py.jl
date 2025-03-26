using PyCall

cd(@__DIR__)

# Import the Python absorb_all function
py"""
import sys
sys.path.append('/Users/C837213770/Desktop/CSU HW/ATS 622/Project 1/src/external scripts')
from absorb import absorb_all as py_absorb_all
"""

# Define the Julia absorb_all function
include("absorb_py.jl")

# Test parameters
frequencies = 0:50:800
temperatures = 250:10:350
pressures = 0.1:50:1000
vapor_pressures = 0:1:15
clw_values = 0:1:15

# Function to compare results
function compare_results(freq, temp, pres, vapor_pres, clw)
    julia_result = absorb_all_unitless(freq, temp, pres, vapor_pres, clw)
    python_result = py"py_absorb_all"(freq, temp, pres, vapor_pres, clw)
    if !isapprox(julia_result, python_result, rtol=1e-10)
        println("Difference found for parameters: freq=$freq, temp=$temp, pres=$pres, vapor_pres=$vapor_pres, clw=$clw")
        println("Julia result: ", julia_result)
        println("Python result: ", python_result)
    end
end

# Iterate over all combinations of parameters
for freq in frequencies
    for temp in temperatures
        for pres in pressures
            for vapor_pres in vapor_pressures
                for clw in clw_values
                    compare_results(freq, temp, pres, vapor_pres, clw)
                end
            end
        end
    end
end

println("Validation completed.")
