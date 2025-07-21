from numpy import sqrt, sin, cos, pi, abs, isnan, zeros


def Interpol(Wp, Wmax, Wmin, Vmax, Vmin):
    """
         ----- Interpolation, RTCA DO-229D (A.4.2.4) -----
         Wp = weighting point
         Wmax = Max weight
         Wmin = Min weight
         Vmax = Value of max weight
         Vmin = Value of min weight

        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021)
    """
    return Vmin + (Vmax - Vmin) * ((Wp - Wmin) / (Wmax - Wmin))


def EstimateTropDalay(ele, lat, h, DOY):
    """
        Estimate tropospheric Delays (A.4.2.4 of RTCA DO-229D)
         Input:
                ele     = Elevation angles in degree
                lat     = latitude of the receiver position on the earth
                h       = Altitude of the receiver position on the earth
                DOY     = Day of year
         Output:
                dTrop   = Tropospheric delays

        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021)
    """
    if isnan(lat) | isnan(h) | isnan(DOY):  # Check input is NaN
        print("Errors: Value of lat, h, or DOY is NaN")

    Lat = abs(lat)

    if lat >= 0:
        Dmin = 28   # for northern latitude
    else:
        Dmin = 211  # for southern latitude

    H = h
    k1 = 77.604     # K/mbar
    k2 = 382000     # K^2/mbar
    Rd = 287.054    # J/(kg.K)
    gm = 9.784      # m/s^2
    g = 9.80665     # m/s^2

    if lat <= 15:
        P0 = 1013.25    # Pressure(mbar)
        T0 = 299.65     # Temperature(K)
        e0 = 26.31      # Water vapor presure(mbar)
        B0 = 6.30e-3    # Temperature lapse rate(K/m)
        Lam0 = 2.77     # Water vapor "lapse rate" (dimensionless)
        Delta_P = 0
        Delta_T = 0
        Delta_e = 0
        Delta_B = 0
        Delta_Lam = 0
    elif 15 < Lat <= 30:
        P0 = Interpol(Lat, 30, 15, 1017.25, 1013.25)
        T0 = Interpol(Lat, 30, 15, 294.15, 299.65)
        e0 = Interpol(Lat, 30, 15, 21.79, 26.31)
        B0 = Interpol(Lat, 30, 15, 6.05e-3, 6.3e-3)
        Lam0 = Interpol(Lat, 30, 15, 3.15, 2.77)
        Delta_P = Interpol(Lat, 30, 15, -3.75, 0)
        Delta_T = Interpol(Lat, 30, 15, 7, 0)
        Delta_e = Interpol(Lat, 30, 15, 8.85, 0)
        Delta_B = Interpol(Lat, 30, 15, 0.25e-3, 0)
        Delta_Lam = Interpol(Lat, 30, 15, 0.33, 0)
    elif 30 < Lat <= 45:
        P0 = Interpol(Lat, 45, 30, 1015.75, 1017.25)
        T0 = Interpol(Lat, 45, 30, 283.15, 294.15)
        e0 = Interpol(Lat, 45, 30, 11.66, 21.79)
        B0 = Interpol(Lat, 45, 30, 5.58e-3, 6.05e-5)
        Lam0 = Interpol(Lat, 45, 30, 2.57, 3.15)
        Delta_P = Interpol(Lat, 45, 30, -2.25, -3.75)
        Delta_T = Interpol(Lat, 45, 30, 11, 7)
        Delta_e = Interpol(Lat, 45, 30, 7.24, 8.85)
        Delta_B = Interpol(Lat, 45, 30, 0.32e-3, 0.25e-3)
        Delta_Lam = Interpol(Lat, 45, 30, 0.46, 0.33)
    elif 45 < Lat <= 60:
        P0 = Interpol(Lat, 60, 45, 1011.75, 1015.75)
        T0 = Interpol(Lat, 60, 45, 272.15, 283.15)
        e0 = Interpol(Lat, 60, 45, 6.78, 11.66)
        B0 = Interpol(Lat, 60, 45, 5.39e-3, 5.58e-3)
        Lam0 = Interpol(Lat, 60, 45, 1.81, 2.57)
        Delta_P = Interpol(Lat, 60, 45, -1.75, -2.25)
        Delta_T = Interpol(Lat, 60, 45, 15, 11)
        Delta_e = Interpol(Lat, 60, 45, 5.36, 7.24)
        Delta_B = Interpol(Lat, 60, 45, 0.81e-3, 0.32e-3)
        Delta_Lam = Interpol(Lat, 60, 45, 0.74, 0.46)
    elif 60 < Lat <= 75:
        P0 = Interpol(Lat, 60, 75, 1013, 1011.75)
        T0 = Interpol(Lat, 60, 75, 263.65, 272.15)
        e0 = Interpol(Lat, 60, 75, 4.11, 6.78)
        B0 = Interpol(Lat, 60, 75, 4.53e-3, 5.39e-3)
        Lam0 = Interpol(Lat, 60, 75, 1.55, 1.81)
        Delta_P = Interpol(Lat, 60, 75, -0.5, -1.75)
        Delta_T = Interpol(Lat, 60, 75, 14.5, 15)
        Delta_e = Interpol(Lat, 60, 75, 3.39, 5.36)
        Delta_B = Interpol(Lat, 60, 75, 0.62e-3, 0.81e-3)
        Delta_Lam = Interpol(Lat, 60, 75, 0.30, 0.74)
    else:   # Lat > 70:
        P0 = 1013
        T0 = 263.65
        e0 = 4.11
        B0 = 4.53e-3
        Lam0 = 1.55
        Delta_P = -0.5
        Delta_T = 14.5
        Delta_e = 3.39
        Delta_B = 0.62e-3
        Delta_Lam = 0.30

    P = P0 - Delta_P * cos(2 * pi * (DOY - Dmin) / 365.25)
    T = T0 - Delta_T * cos(2 * pi * (DOY - Dmin) / 365.25)
    e = e0 - Delta_e * cos(2 * pi * (DOY - Dmin) / 365.25)
    B = B0 - Delta_B * cos(2 * pi * (DOY - Dmin) / 365.25)
    Lam = Lam0 - Delta_Lam * cos(2 * pi * (DOY - Dmin) / 365.25)
    Zhyd = ((10 ** -6) * k1 * Rd * P) / gm
    Zwet = (((10 ** -6) * k2 * Rd) / (gm * (Lam + 1) - B0 * Rd)) * (e / T)
    d_hyd = ((1 - (B * H / T)) ** (g / (Rd * B))) * Zhyd    # Range delays for the hydrostatic equilibrium
    d_wet = ((1 - (B * H / T)) ** (((Lam + 1) * g / (Rd * B)) - 1)) * Zwet  # Range delays for the water vapor

    return (d_hyd + d_wet) * (1.001 / sqrt(0.002001 + sin(ele * pi / 180) ** 2))


def Klobuchar(alpha, beta, ele, azi, lat, lon, t_local):
    """
        Estimate the ionospheric delay based on the Klobuchar model (20.3.3.5.2.5 of IS-GPS-200D).
         Input:
                alpha   = ION ALPHA
                beta    = ION BETA
                ele     = Elevation angle in degree
                azi     = Azimuth angle in degree
                lat     = latitude in degree
                lon     = Longitude in degree
                SOD     = Second of day
        Output:
                Ionospheric delay estimations

        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021).
    """

    # --- Ionospheric model (20.3.3.5.2.5 of IS-GPS-200D) ---
    El = ele / 180      # Elevation angle (semi-circles)
    Az = azi / 180      # Azimuth angle (semi-circles)
    Shi = 0.0137 / (El + 0.11) - 0.022      # Earth's central angle
    Lam_u = lon / 180       # User geodetic longitude(semi-circles) WGS-84
    Phi_u = lat / 180       # User geodetic latitude(semi-circles) WGS-84
    Phi_i = Phi_u + Shi * cos(Az * pi)      # Geodetic latitude (semi-circles)
    Phi_i = (Phi_i > 0.416) * 0.416 + (Phi_i < -0.416) * (-0.416) + (abs(Phi_i) <= 0.416) * Phi_i
    Lam_i = Lam_u + Shi * (sin(Az * pi) / cos(Phi_i * pi))      # Geodetic longitude (semi-circles)
    Phi_m = Phi_i + 0.064 * cos((Lam_i - 1.617) * pi)       # Geomagnetic latitude(assumed height 350 km) (semi-circles)
    t_m = (Lam_i * 4.32e4) + t_local      # Local time(sec)
    t_m = (t_m >= 86400) * (t_m - 86400) + (t_m < 0) * (t_m + 86400) + (t_m < 86400) * (t_m >= 0) * t_m
    F = 1 + 16 * (0.53 - El) ** 3   # Obliquity factor
    PER = zeros((len(Phi_m), 1))
    AMP = zeros((len(Phi_m), 1))

    for i in range(len(Phi_m)):
        PER[i] = beta[0] + beta[1] * Phi_m[i] + beta[2] * Phi_m[i] ** 2 + beta[3] * Phi_m[i] ** 3
        AMP[i] = alpha[0] + alpha[1] * Phi_m[i] + alpha[2] * Phi_m[i] ** 2 + alpha[3] * Phi_m[i] ** 3

    PER = (PER < 72000) * 72000 + (PER >= 72000) * PER
    AMP = (AMP >= 0) * AMP
    Xm = 2 * pi * (t_m - 50400) / PER
    return (abs(Xm) <= 1.57) * (F * (5e-9 + AMP * (1 - (Xm ** 2) / 2 + (Xm ** 4)/24))) + (abs(Xm) > 1.57) * (F * 5e-9)
