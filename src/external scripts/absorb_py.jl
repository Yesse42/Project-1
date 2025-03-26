using LinearAlgebra, Unitful

function absorb_all(freq::Unitful.Frequency, temp::Unitful.Temperature, pres::Unitful.Pressure, vapor_pres::Unitful.Pressure, clw::Unitful.Density)
    """
    Wrapper function to compute the absorption coefficient for water vapor (H2O),
    oxygen (O2), and nitrogen (N2) with unit handling.

    Parameters:
    - freq: Frequency [GHz]
    - temp: Layer average temperature [K]
    - pres: Layer average pressure [mbar]
    - vapor_pres: Layer average water vapor pressure [mbar]
    - clw: Cloud liquid water density [g/m^3]

    Returns:
    - kabs_clear: Absorption coefficient [km^-1]
    """
    # Convert frequency to GHz and strip units
    freqy_val = ustrip(uconvert(u"GHz", freq))

    # Check frequency is in the 0-800GHz range
    if freqy_val < 0 || freqy_val > 800
        error("Frequency must be in the range 0-800 GHz")
    end
    
    # Convert to specified units and strip units
    temp_val = ustrip(uconvert(u"K", temp))
    pres_val = ustrip(uconvert(u"mbar", pres))
    vapor_pres_val = ustrip(uconvert(u"mbar", vapor_pres))
    clw_val = ustrip(uconvert(u"g/m^3", clw))

    # Compute the absorption coefficient without units
    kabs_all_unitless = absorb_all_unitless(freqy_val, temp_val, pres_val, vapor_pres_val, clw_val)

    return kabs_all_unitless * u"km^-1"
end

function absorb_all_unitless(freqy, temp, pres, vapor_pres, clw)
    """
    Function to compute the absorption coefficient for water vapor (H2O),
    oxygen (O2), and nitrogen (N2) without unit handling.

    Parameters:
    - freqy: Frequency [GHz]
    - temp: Layer average temperature [K]
    - pres: Layer average pressure [mb]
    - vapor_pres: Layer average water vapor pressure [mb]
    - clw: Cloud liquid water density [g/m^3]

    Returns:
    - kabs_clear: Absorption coefficient [km^-1]
    """
    # Compute vapor density, rho [g/m^3]
    rho = vapor_pres * 100 * 18 / (8.314 * temp)

    # Calculate individual absorption coefficients
    kabs_H2O = abs_H2O(temp, pres, rho, freqy)
    kabs_O2 = abs_O2(temp, pres, rho, freqy)
    kabs_N2 = abs_N2(temp, pres, freqy)
    kabs_clw = absorb_clw(freqy, temp, clw)

    # Total absorption coefficient
    kabs_all = kabs_H2O + kabs_O2 + kabs_N2 + kabs_clw

    return kabs_all
end

function abs_H2O(T, P, RHO, F)
    """
    Compute absorption coefficient in the atmosphere due to water vapor.

    Parameters:
    T : Float64
        Temperature in Kelvin
    P : Float64
        Pressure in millibar (0.1 to 1000)
    RHO : Float64
        Water vapor density in g/m^3
    F : Float64
        Frequency in GHz (0 to 800)

    Returns:
    ABH2O : Float64
        Absorption coefficient in Nepers/Km
    """
    FL = [22.2351, 183.3101, 321.2256, 325.1529, 380.1974, 439.1508, 443.0183, 448.0011, 470.8890, 474.6891, 488.4911, 556.9360, 620.7008, 752.0332, 916.1712]
    S1 = [0.1314E-13, 0.2279E-11, 0.8058E-13, 0.2701E-11, 0.2444E-10, 0.2185E-11, 0.4637E-12, 0.2568E-10, 0.8392E-12, 0.3272E-11, 0.6676E-12, 0.1535E-08, 0.1711E-10, 0.1014E-08, 0.4238E-10]
    B2 = [2.144, 0.668, 6.179, 1.541, 1.048, 3.595, 5.048, 1.405, 3.597, 2.379, 2.852, 0.159, 2.391, 0.396, 1.441]
    W3 = [0.00281, 0.00287, 0.0023, 0.00278, 0.00287, 0.0021, 0.00186, 0.00263, 0.00215, 0.00236, 0.0026, 0.00321, 0.00244, 0.00306, 0.00267]
    X = [0.69, 0.64, 0.67, 0.68, 0.54, 0.63, 0.60, 0.66, 0.66, 0.65, 0.69, 0.69, 0.71, 0.68, 0.70]
    WS = [0.01349, 0.01491, 0.0108, 0.0135, 0.01541, 0.0090, 0.00788, 0.01275, 0.00983, 0.01095, 0.01313, 0.01320, 0.01140, 0.01253, 0.01275]
    XS = [0.61, 0.85, 0.54, 0.74, 0.89, 0.52, 0.50, 0.67, 0.65, 0.64, 0.72, 1.0, 0.68, 0.84, 0.78]
    SR = [0.0, -0.017, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    DF = [0.0, 0.0]

    if RHO <= 0
        return 0.0
    end

    PVAP = RHO * T / 217.0
    PDA = P - PVAP
    DEN = 3.335E16 * RHO
    TI = 300.0 / T
    TI2 = TI ^ 2.5

    CON = (5.43E-10 * PDA * TI^3 + 1.8E-8 * PVAP * TI^7.5) * PVAP * F^2

    SUM = 0.0
    for I in 1:15
        WIDTH = W3[I] * PDA * TI ^ X[I] + WS[I] * PVAP * TI ^ XS[I]
        SHIFT = SR[I] * WIDTH
        WSQ = WIDTH ^ 2
        S = S1[I] * TI2 * exp(B2[I] * (1.0 - TI))

        DF[1] = F - FL[I] - SHIFT
        DF[2] = F + FL[I] + SHIFT

        BASE = WIDTH / (562500.0 + WSQ)

        RES = 0.0
        for J in 1:2
            if abs(DF[J]) < 750.0
                RES += WIDTH / (DF[J]^2 + WSQ) - BASE
            end
        end

        SUM += S * RES * (F / FL[I]) ^ 2
    end

    ABH2O = 0.3183E-4 * DEN * SUM + CON
    return ABH2O
end

function abs_O2(TEMP, PRES, VAPDEN, FREQY)
    F = [118.7503, 56.2648, 62.4863, 58.4466, 60.3061, 59.5910, 59.1642, 60.4348, 58.3239, 61.1506, 57.6125, 61.8002, 56.9682, 62.4112, 56.3634, 62.9980, 55.7838, 63.5685, 55.2214, 64.1278, 54.6712, 64.6789, 54.1300, 65.2241, 53.5957, 65.7648, 53.0669, 66.3021, 52.5424, 66.8368, 52.0214, 67.3696, 51.5034, 67.9009, 368.4984, 424.7632, 487.2494, 715.3931, 773.8397, 834.1458]
    S300 = [0.2936E-14, 0.8079E-15, 0.2480E-14, 0.2228E-14, 0.3351E-14, 0.3292E-14, 0.3721E-14, 0.3891E-14, 0.3640E-14, 0.4005E-14, 0.3227E-14, 0.3715E-14, 0.2627E-14, 0.3156E-14, 0.1982E-14, 0.2477E-14, 0.1391E-14, 0.1808E-14, 0.9124E-15, 0.1230E-14, 0.5603E-15, 0.7842E-15, 0.3228E-15, 0.4689E-15, 0.1748E-15, 0.2632E-15, 0.8898E-16, 0.1389E-15, 0.4264E-16, 0.6899E-16, 0.1924E-16, 0.3229E-16, 0.8191E-17, 0.1423E-16, 0.6494E-15, 0.7083E-14, 0.3025E-14, 0.1835E-14, 0.1158E-13, 0.3993E-14]
    BE = [0.009, 0.015, 0.083, 0.084, 0.212, 0.212, 0.391, 0.391, 0.626, 0.626, 0.915, 0.915, 1.260, 1.260, 1.660, 1.665, 2.119, 2.115, 2.624, 2.625, 3.194, 3.194, 3.814, 3.814, 4.484, 4.484, 5.224, 5.224, 6.004, 6.004, 6.844, 6.844, 7.744, 7.744, 0.048, 0.044, 0.049, 0.145, 0.141, 0.145]
    WB300 = 0.56
    X = 0.8
    W300 = [1.67, 1.646, 1.468, 1.449, 1.382, 1.360, 1.319, 1.297, 1.266, 1.248, 1.221, 1.207, 1.181, 1.171, 1.144, 1.139, 1.110, 1.108, 1.079, 1.078, 1.05, 1.05, 1.02, 1.02, 1.00, 1.00, 0.97, 0.97, 0.94, 0.94, 0.92, 0.92, 0.89, 0.89, 1.64, 1.64, 1.64, 1.81, 1.81, 1.81]
    Y300 = [-0.036, 0.2408, -0.3486, 0.5227, -0.5430, 0.5877, -0.3970, 0.3237, -0.1348, 0.0311, 0.0725, -0.1663, 0.2832, -0.3629, 0.3970, -0.4599, 0.4695, -0.5199, 0.5187, -0.5597, 0.5903, -0.6246, 0.6656, -0.6942, 0.7086, -0.7325, 0.7348, -0.7546, 0.7702, -0.7864, 0.8083, -0.8210, 0.8439, -0.8529, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    V = [0.0079, -0.0978, 0.0844, -0.1273, 0.0699, -0.0776, 0.2309, -0.2825, 0.0436, -0.0584, 0.6056, -0.6619, 0.6451, -0.6759, 0.6547, -0.6675, 0.6135, -0.6139, 0.2952, -0.2895, 0.2654, -0.2590, 0.3750, -0.3680, 0.5085, -0.5002, 0.6206, -0.6091, 0.6526, -0.6393, 0.6640, -0.6475, 0.6729, -0.6545, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    TH = 300.0 / TEMP
    TH1 = TH - 1.0
    B = TH ^ X
    PRESWV = VAPDEN * TEMP / 217.0
    PRESDA = PRES - PRESWV
    DEN = 0.001 * (PRESDA * B + 1.1 * PRESWV * TH)
    DENS = 0.001 * (PRESDA * TH ^ 0.9 + 1.1 * PRESWV * TH)
    DFNR = WB300 * DEN
    SUM = 1.6e-17 * FREQY ^ 2 * DFNR / (TH * (FREQY ^ 2 + DFNR ^ 2))

    for K in 1:40
        if K == 1
            DF = W300[1] * DENS
        else
            DF = W300[K] * DEN
        end

        FCEN = F[K]
        Y = 0.001 * PRES * B * (Y300[K] + V[K] * TH1)
        STR = S300[K] * exp(-BE[K] * TH1)
        SF1 = (DF + (FREQY - FCEN) * Y) / ((FREQY - FCEN) ^ 2 + DF ^ 2)
        SF2 = (DF - (FREQY + FCEN) * Y) / ((FREQY + FCEN) ^ 2 + DF ^ 2)
        SUM += STR * (SF1 + SF2) * (FREQY / F[K]) ^ 2
    end

    O2ABS = 0.5034e12 * SUM * PRESDA * TH ^ 3 / π
    O2ABS = max(O2ABS, 0.0)
    return O2ABS
end

function abs_N2(temp, pres, freqy)
    """
    Function to calculate the nitrogen absorption coefficient (N2).

    Parameters:
    - temp: Temperature [K]
    - pres: Pressure [mb]
    - freqy: Frequency [GHz]

    Returns:
    - kabs_N2: Absorption coefficient for N2 [km^-1]
    """
    th = 300.0 / temp
    fdep = 0.5 + 0.5 / (1.0 + (freqy / 450.0) ^ 2)
    bf = 6.5E-14 * fdep * pres ^ 2 * freqy ^ 2 * th ^ 3.6
    absn2 = 1.29 * bf
    return absn2
end

function absorb_clw(freq, temp, water)
    """
    Computes absorption in nepers/km by suspended water droplets

    Arguments (Input):
    - water: Water in g/m^3
    - freq: Frequency in GHz (valid from 0 to 1000 GHz)
    - temp: Temperature in Kelvin

    Returns:
    - kabs_clw: Absorption coefficient for cloud liquid water [km^-1]
    """
    if water <= 0
        return 0.0
    end

    THETA1 = 1.0 - 300.0 / temp
    EPS0 = 77.66 - 103.3 * THETA1
    EPS1 = 0.0671 * EPS0
    EPS2 = 3.52
    FP = 20.1 * exp(7.88 * THETA1)
    FS = 39.8 * FP
    EPS = (EPS0 - EPS1) / complex(1.0, freq / FP) + (EPS1 - EPS2) / complex(1.0, freq / FS) + EPS2
    RE = (EPS - 1.0) / (EPS + 2.0)
    kabs_clw = -0.06286 * imag(RE) * freq * water

    return kabs_clw
end