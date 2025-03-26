using DataFrames, Unitful

data = [
    (0.0u"km", 1.013e3u"mbar", 288.1u"K", 1.225e3u"g/m^3", 5.9e0u"g/m^3", 5.4e-5u"g/m^3"),
    (1.0u"km", 8.986e2u"mbar", 281.6u"K", 1.111e3u"g/m^3", 4.2e0u"g/m^3", 5.4e-5u"g/m^3"),
    (2.0u"km", 7.950e2u"mbar", 275.1u"K", 1.007e3u"g/m^3", 2.9e0u"g/m^3", 5.4e-5u"g/m^3"),
    (3.0u"km", 7.012e2u"mbar", 268.7u"K", 9.093e2u"g/m^3", 1.8e0u"g/m^3", 5.0e-5u"g/m^3"),
    (4.0u"km", 6.166e2u"mbar", 262.2u"K", 8.193e2u"g/m^3", 1.1e0u"g/m^3", 4.6e-5u"g/m^3"),
    (5.0u"km", 5.405e2u"mbar", 255.7u"K", 7.364e2u"g/m^3", 6.4e-1u"g/m^3", 4.5e-5u"g/m^3"),
    (6.0u"km", 4.722e2u"mbar", 249.2u"K", 6.601e2u"g/m^3", 3.8e-1u"g/m^3", 4.5e-5u"g/m^3"),
    (7.0u"km", 4.111e2u"mbar", 242.7u"K", 5.900e2u"g/m^3", 2.1e-1u"g/m^3", 4.8e-5u"g/m^3"),
    (8.0u"km", 3.565e2u"mbar", 236.2u"K", 5.258e2u"g/m^3", 1.2e-1u"g/m^3", 5.2e-5u"g/m^3"),
    (9.0u"km", 3.080e2u"mbar", 229.7u"K", 4.671e2u"g/m^3", 4.6e-2u"g/m^3", 7.1e-5u"g/m^3"),
    (10.0u"km", 2.650e2u"mbar", 223.2u"K", 4.135e2u"g/m^3", 1.8e-2u"g/m^3", 9.0e-5u"g/m^3"),
    (11.0u"km", 2.270e2u"mbar", 216.8u"K", 3.648e2u"g/m^3", 8.2e-3u"g/m^3", 1.3e-4u"g/m^3"),
    (12.0u"km", 1.940e2u"mbar", 216.6u"K", 3.119e2u"g/m^3", 3.7e-3u"g/m^3", 1.6e-4u"g/m^3"),
    (13.0u"km", 1.658e2u"mbar", 216.6u"K", 2.666e2u"g/m^3", 1.8e-3u"g/m^3", 1.7e-4u"g/m^3"),
    (14.0u"km", 1.417e2u"mbar", 216.6u"K", 2.279e2u"g/m^3", 8.4e-4u"g/m^3", 1.9e-4u"g/m^3"),
    (15.0u"km", 1.211e2u"mbar", 216.6u"K", 1.948e2u"g/m^3", 7.2e-4u"g/m^3", 2.1e-4u"g/m^3"),
    (16.0u"km", 1.035e2u"mbar", 216.6u"K", 1.665e2u"g/m^3", 6.1e-4u"g/m^3", 2.3e-4u"g/m^3"),
    (17.0u"km", 8.850e1u"mbar", 216.6u"K", 1.423e2u"g/m^3", 5.2e-4u"g/m^3", 2.8e-4u"g/m^3"),
    (18.0u"km", 7.565e1u"mbar", 216.6u"K", 1.216e2u"g/m^3", 4.4e-4u"g/m^3", 3.2e-4u"g/m^3"),
    (19.0u"km", 6.467e1u"mbar", 216.6u"K", 1.040e2u"g/m^3", 4.4e-4u"g/m^3", 3.5e-4u"g/m^3"),
    (20.0u"km", 5.529e1u"mbar", 216.6u"K", 8.891e1u"g/m^3", 4.4e-4u"g/m^3", 3.8e-4u"g/m^3"),
    (21.0u"km", 4.729e1u"mbar", 217.6u"K", 7.572e1u"g/m^3", 4.8e-4u"g/m^3", 3.8e-4u"g/m^3"),
    (22.0u"km", 4.047e1u"mbar", 218.6u"K", 6.451e1u"g/m^3", 5.2e-4u"g/m^3", 3.9e-4u"g/m^3"),
    (23.0u"km", 3.467e1u"mbar", 218.6u"K", 5.500e1u"g/m^3", 5.7e-4u"g/m^3", 3.8e-4u"g/m^3"),
    (24.0u"km", 2.972e1u"mbar", 220.6u"K", 4.694e1u"g/m^3", 6.1e-4u"g/m^3", 3.6e-4u"g/m^3"),
    (25.0u"km", 2.549e1u"mbar", 221.6u"K", 4.008e1u"g/m^3", 6.6e-4u"g/m^3", 3.4e-4u"g/m^3"),
    (30.0u"km", 1.197e1u"mbar", 226.5u"K", 1.841e1u"g/m^3", 3.8e-4u"g/m^3", 2.0e-4u"g/m^3"),
    (35.0u"km", 5.746e0u"mbar", 236.5u"K", 8.463e0u"g/m^3", 1.6e-4u"g/m^3", 1.1e-4u"g/m^3"),
    (40.0u"km", 2.871e0u"mbar", 250.4u"K", 3.996e0u"g/m^3", 6.7e-5u"g/m^3", 4.9e-5u"g/m^3"),
    (45.0u"km", 1.491e0u"mbar", 264.2u"K", 1.996e0u"g/m^3", 3.2e-5u"g/m^3", 1.7e-5u"g/m^3"),
    (50.0u"km", 7.978e-1u"mbar", 270.6u"K", 1.027e0u"g/m^3", 1.2e-5u"g/m^3", 4.0e-6u"g/m^3"),
    (75.0u"km", 5.520e-2u"mbar", 219.7u"K", 8.754e-2u"g/m^3", 1.5e-7u"g/m^3", 8.6e-8u"g/m^3"),
    (100.0u"km", 3.008e-4u"mbar", 210.0u"K", 4.989e-4u"g/m^3", 1.0e-9u"g/m^3", 4.3e-11u"g/m^3")
]

standard_atm = DataFrame(data, [:Hgt, :Pressure, :T, :density, :H20, :O3])

#convert the water vapor density to a vapor Pressure
standard_atm.H20 = uconvert.(u"hPa", standard_atm.H20 ./ 18.01528u"g/mol" .* 8.31446261815324u"J/(mol*K)" .* standard_atm.T)

standard_atm;

