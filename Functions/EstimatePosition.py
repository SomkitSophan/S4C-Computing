from numpy import array, tan, arctan2, arctan, sqrt, sin, arcsin, cos, linalg, pi, where, deg2rad, rad2deg
from numpy import isnan, zeros, ones, dot, sum, append, inf
from Functions.AtmosphereDelays import EstimateTropDalay, Klobuchar
from Functions.Coordinate import ecef2enu, ecef2lla, Rotation_Z
import datetime
import operator

def GPS_xyz(nav, epoch, SOD, R_raw, PRN, XYZ_r):
    """"
        Calculate GPS position
         Inputs:
                nav     = Ephemeris data
                epoch   = Date of observations
                SOD     = Second of day
                R_raw   = Pseudorange measurement
                PRN     = PRN number
                R_xyz   = Receiver position (ECEF) (m)
         Outputs:
               XYZ_s    = Satellite position (m)
               t_biasL1 = Satellite clock bias (sec)
               ele      = Elevation angle (deg)
               azi      = Azimuth angle (deg)
        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021)
    """
    # Constant
    gm = 3.9860050 * 10 ** 14  # Earth's universal gravitational parameter (m^3/s^2)
    We = 7.2921151467 * 10 ** -5  # Earth rotation rate (rad/sec)
    c = 299792458  # Speed of light (m/s)

    y = epoch.year  # Year
    m = epoch.month  # Month
    d = epoch.day  # Day of month

    if datetime.date(y, m, d).isoweekday() < 7:
        SOW = (datetime.date(y, m, d).isoweekday()) * (24 * 60 * 60)  # Computed the second of the GPS Week (SOW)
    else:
        SOW = 0

    # --- Estimated the time to compute the satellite position ---
    Tr = R_raw / c
    Ts = SOW + SOD - Tr

    # --- Read Ephemeris data & Find the accurate ephemeris ---
    ind = isnan((nav['IODE'][PRN])) == False
    if ind.size == 1:
        # Orbit Parameters
        a = nav['sqrtA'][PRN] ** 2  # Semi-major axis  (m)
        e = nav['e'][PRN]  # Eccentricity
        w0 = nav['omega'][PRN]  # Argument of perigee (rad)
        W0 = nav['OMEGA0'][PRN]  # Right ascension of ascending node (rad)
        Wdot = nav['OMEGAdot'][PRN]  # Rate of right a
        i0 = nav['I0'][PRN]  # Inclination (rad)
        idot = nav['Idot'][PRN]  # Rate of inclination (rad/sec)
        M0 = nav['M0'][PRN] # Mean anomaly (rad)
        delta_n = nav['delta_n'][PRN]  # Mean motion rate (rad/sec)

        # Correction coefficients
        Cuc = nav['cuc'][PRN]  # Argument of perigee (cos)   (rad)
        Cus = nav['cus'][PRN]  # Argument of perigee (sine)  (rad)
        Crc = nav['crc'][PRN]  # Orbit radius        (cos)   (m)
        Crs = nav['crs'][PRN]  # Orbit radius        (sine)  (m)
        Cic = nav['cic'][PRN]  # Inclination         (cos)   (rad)
        Cis = nav['cis'][PRN]  # Inclination         (sine)  (rad)

        # Time
        Toe = nav['Toe'][PRN]  # Time of Ephemeris (SOW : sec of GPS week)
        GPS_week = nav['SOW'][PRN]  # GPS Week
        Ttm = nav['Ts'][PRN]  # Transmission time of message-604800  (SOW : sec of GPS week)

        # Clock
        af0 = nav['af0'][PRN]  # Clock Bias (sec)
        af1 = nav['af1'][PRN]  # Clock Drift (sec/sec)
        af2 = nav['af2'][PRN]  # Clock Drift rate (sec/sec^2)
        Tgd = nav['TGD'][PRN]  # Time Group delay (sec)

        # Status
        SV_health = nav['SVh'][PRN]  # SV Health
        SV_accuracy = nav['SVa'][PRN]  # SV Accuracy
        L2_P_flag = nav['flagP2'][PRN]  # L2 P data flag
        L2_code = nav['code2'][PRN]  # Code on L2 channel
        IODC = nav['IODC'][PRN] # Issue of Data, Clock
        IODE = nav['IODE'][PRN]  # Issue of Data, Ephemeris
    else:
        Toe = nav['Toe'][PRN].loc[ind].values  # Time of Ephemeris (SOW : sec of GPS week)
        min_ind, min_value = min(enumerate(abs(Toe - Ts)), key=operator.itemgetter(1))

        # Orbit Parameters
        a = nav['sqrtA'][PRN].loc[ind].values[min_ind] ** 2  # Semi-major axis  (m)
        e = nav['e'][PRN].loc[ind].values[min_ind]  # Eccentricity
        w0 = nav['omega'][PRN].loc[ind].values[min_ind]  # Argument of perigee (rad)
        W0 = nav['OMEGA0'][PRN].loc[ind].values[min_ind]  # Right ascension of ascending node (rad)
        Wdot = nav['OMEGAdot'][PRN].loc[ind].values[min_ind]  # Rate of right a
        i0 = nav['I0'][PRN].loc[ind].values[min_ind]  # Inclination (rad)
        idot = nav['Idot'][PRN].loc[ind].values[min_ind]  # Rate of inclination (rad/sec)
        M0 = nav['M0'][PRN].loc[ind].values[min_ind]  # Mean anomaly (rad)
        delta_n = nav['delta_n'][PRN].loc[ind].values[min_ind]  # Mean motion rate (rad/sec)

        # Correction coefficients
        Cuc = nav['cuc'][PRN].loc[ind].values[min_ind]  # Argument of perigee (cos)   (rad)
        Cus = nav['cus'][PRN].loc[ind].values[min_ind]  # Argument of perigee (sine)  (rad)
        Crc = nav['crc'][PRN].loc[ind].values[min_ind]  # Orbit radius        (cos)   (m)
        Crs = nav['crs'][PRN].loc[ind].values[min_ind]  # Orbit radius        (sine)  (m)
        Cic = nav['cic'][PRN].loc[ind].values[min_ind]  # Inclination         (cos)   (rad)
        Cis = nav['cis'][PRN].loc[ind].values[min_ind]  # Inclination         (sine)  (rad)

        # Time
        Toe = nav['Toe'][PRN].loc[ind].values[min_ind]  # Time of Ephemeris (SOW : sec of GPS week)
        GPS_week = nav['SOW'][PRN].loc[ind].values[min_ind]  # GPS Week
        Ttm = nav['Ts'][PRN].loc[ind].values[min_ind]  # Transmission time of message-604800  (SOW : sec of GPS week)

        # Clock
        af0 = nav['af0'][PRN].loc[ind].values[min_ind]  # Clock Bias (sec)
        af1 = nav['af1'][PRN].loc[ind].values[min_ind]  # Clock Drift (sec/sec)
        af2 = nav['af2'][PRN].loc[ind].values[min_ind]  # Clock Drift rate (sec/sec^2)
        Tgd = nav['TGD'][PRN].loc[ind].values[min_ind]  # Time Group delay (sec)

        # Status
        SV_health = nav['SVh'][PRN].loc[ind].values[min_ind]  # SV Health
        SV_accuracy = nav['SVa'][PRN].loc[ind].values[min_ind]  # SV Accuracy
        L2_P_flag = nav['flagP2'][PRN].loc[ind].values[min_ind]  # L2 P data flag
        L2_code = nav['code2'][PRN].loc[ind].values[min_ind]  # Code on L2 channel
        IODC = nav['IODC'][PRN].loc[ind].values[min_ind]  # Issue of Data, Clock
        IODE = nav['IODE'][PRN].loc[ind].values[min_ind]  # Issue of Data, Ephemeris

    Tk = Ts - Toe  # Time elapsed since Toe(SOW : sec of GPS week)
    M = M0 + (sqrt(gm / a ** 3) + delta_n) * Tk  # Mean anomaly at Tk

    # Iterative solution for E
    E_old = M
    dE = 1

    while dE > 10 ** -12:       # A predicted loop of the Eccentric anomaly
        E = E_old - (E_old - e * sin(E_old) - M) / (1 - e * cos(E_old))  # Eccentric anomaly
        dE = abs(E - E_old)
        E_old = E

    v = arctan2(sqrt(1 - e ** 2) * sin(E), cos(E) - e)  # True anomaly
    p = v + w0  # Argument of perigee + True anomaly or Argument of latitude

    # Corrected the orbital errors
    u = p + Cuc * cos(2 * p) + Cus * sin(2 * p)  # Corrected argument of latitude(rad)
    r = a * (1 - e * cos(E)) + Crc * cos(2 * p) + Crs * sin(2 * p)  # Corrected radius(meter)
    i = i0 + idot * Tk + Cic * cos(2 * p) + Cis * sin(2 * p)  # Corrected inclination(rad)

    W = W0 + (Wdot - We) * Tk - (We * Toe)  # Argument of ascending node

    # rotation matrix
    R = array([cos(W) * cos(u) - sin(W) * sin(u) * cos(i),
               sin(W) * cos(u) + cos(W) * sin(u) * cos(i),
               sin(u) * sin(i)])

    XYZ_s = r * R  # Satellite position vector (ECEF with WGS-84)

    # --- Computed the clock errors ---
    F = - 4.442807633e-10  # Constant  (sec/meter^1/2)
    Delta_tr = F * e * sqrt(a) * sin(E)  # Relativistic correction term
    Delta_tsv = af0 + af1 * Tk + af2 * Tk ** 2 + Delta_tr - Tgd  # Delta clock
    dDelta_tsv = af1 + 2 * af2 * Tk  # Clock drift
    t_biasL1 = Delta_tsv + dDelta_tsv  # The satellite clock bias of the L1 frequency

    # --- The ENU = > Elevation & Azimuth ---
    [east, north, up] = ecef2enu(XYZ_r, XYZ_s)
    ele = arcsin(up) * 180 / pi  # Elevation angle(degree)
    azi = arctan2(east, north) * 180 / pi  # Azimuth angle (degree)

    return XYZ_s.reshape(1, -1), t_biasL1, ele, azi


def EstimateReceiverPosition(nav, epoch, Rn, PRN, ele_mask, XYZ_r, sec, K):
    """"
            Calculate Receiver position
             Inputs:
                    nav         = Ephemeris data
                    epoch       = Date of observations
                    Rn          = Pseudorange
                    PRN         = PRN number
                    ele_mask    = Elevation cut-off
                    XYZ_r       = Receiver position (ECEF) (m)
                    sec         = Sampling time
                    K           = Iterative count
             Outputs:
                   XYZ_r        = Receiver position (m)
                   LLA_r        = Latitude & Longitude (degree), Altitude (m)
                   Bias         = Bias (m)
            CSSRG-LAB of KMITL, Thailand.
            Version 1 by Somkit Sophan (July 2021)
        """

    c = 299792458  # Speed of light (m/s)
    Nav = nav[0].ephemeris  # Read Navigation data
    Ionprm = nav[0].Ion  # Coefficients of Iono-delays (Klobuchar model)
    XYZ_r = XYZ_r.reshape(-1, 1)
    LLA_u = ecef2lla(XYZ_r)
    DOY = epoch.timetuple().tm_yday  # Day-of-Year
    Flag = 1
    A = []

    PRN_len = len(PRN)
    if PRN_len < 4:  # Check to need at least 4 satellites
        print("GPS data less than 4 satellites")
        XYZ_r = zeros((3, 1))
        LLA_r = zeros((3, 1))
        Bias = 0
        HDOP = 0
        VDOP = 0
        Flag = 0
    else:
        SOD = sec + epoch.hour * 60 * 60 + epoch.minute * 60 + epoch.second  # Second-of-Day
        XYZ_s = zeros((PRN_len, 3))  # Satellite position variable
        Clock_s = zeros((PRN_len, 1))  # Satellite-clock offset variable
        Ele = zeros((PRN_len, 1))  # Elevation angle variable
        Azi = zeros((PRN_len, 1))  # Azimuth angle variable
        t_gps = zeros((PRN_len, 1))  # GPS time

        for ind in range(PRN_len):  # A predicted loop of the satellite positions
            if not isnan(Rn[ind]):
                for i in range(0, nav.__len__(), 1):
                    Buff = nav[i].ephemeris
                    try:
                        if sum(isnan(Buff['IODE'][PRN[ind]]) == False) > 0:
                            Nav = nav[i].ephemeris  # Read Navigation data
                        break
                    except:
                        continue

                try:
                    XYZ_s[ind], Clock_s[ind], Ele[ind], Azi[ind] = GPS_xyz(Nav, epoch, SOD, Rn[ind], PRN[ind], XYZ_r)
                    t_gps[ind] = SOD - Rn[ind] / c
                except:
                    continue

            else:
                continue

        # Filter by zero values
        Ind_cut = where((t_gps == 0) == False)
        Rn = Rn[Ind_cut[0]]
        XYZ_s = XYZ_s[Ind_cut[0]]
        Clock_s = Clock_s[Ind_cut[0]]
        Ele = Ele[Ind_cut[0]]
        Azi = Azi[Ind_cut[0]]
        t_gps = t_gps[Ind_cut[0]]
        PRN = PRN[Ind_cut[0]]
        PRN_len = len(PRN)

        if K > 0:
            # Filter by ele_mask
            Ind_cut = where((Ele < ele_mask | isnan(Ele) | (Clock_s == 0)) == False)
            Rn = Rn[Ind_cut[0]]
            XYZ_s = XYZ_s[Ind_cut[0]]
            Clock_s = Clock_s[Ind_cut[0]]
            Ele = Ele[Ind_cut[0]]
            Azi = Azi[Ind_cut[0]]
            t_gps = t_gps[Ind_cut[0]]
            PRN = PRN[Ind_cut[0]]
            PRN_len = len(PRN)

            # Iono-delay predictions
            dIon = Klobuchar(Ionprm[0], Ionprm[1], Ele, Azi, LLA_u[0], LLA_u[1], t_gps)
            IC = c * dIon

            # Tropo-delay predictions
            TC = EstimateTropDalay(Ele, LLA_u[0], LLA_u[2], DOY)

        else:
            IC = zeros((PRN_len, 1))
            TC = zeros((PRN_len, 1))

        if (PRN_len < 4) | (sum(XYZ_s) < 1):  # Check to need at least 4 satellites
            print("GPS data less than 4 satellites")
            XYZ_r = zeros((3, 1))
            LLA_r = zeros((3, 1))
            Bias = 0
            HDOP = 0
            VDOP = 0
            Flag = 0
        else:
            # Pseudorange estimations
            R = Rn + (c * Clock_s) - TC - IC

            # Sagnac effect corrections
            XYZ_s_new = zeros((PRN_len, 3))  # New satellite position variable
            for ind in range(PRN_len):  # A corrected loop of sagnac effects
                t_Rota = R[ind] / c
                XYZ_s_new[ind] = Rotation_Z(XYZ_s[ind], t_Rota)

            begin_K = 0  # Iterative count
            end_K = 100  # Iterative number
            Bias = 0  # initial bias
            dX = inf  # Initial dX

            while (linalg.norm(dX) > 10 ** -3) and begin_K < end_K:  # Iterative LSE loop
                begin_K += 1
                # Calculated the line-of-sight vectors and ranges from satellite to receiver
                dXYZ = XYZ_s_new - dot(ones((PRN_len, 1)), XYZ_r.reshape(1, -1))
                R_dXYZ = sqrt(sum(dXYZ ** 2, axis=1)).reshape(-1, 1)
                abc = - dXYZ / dot(R_dXYZ, ones((1, 3)))  # Line of sight unit vectors

                # Calculated the a-priori range
                R_hat = R_dXYZ - Bias
                Z = R - R_hat  # Matrix Z
                A = append(abc, - ones((PRN_len, 1)), axis=1)  # Matrix A

                # LSE solution
                dX = dot(dot(linalg.inv(dot(A.T, A)), A.T), Z)

                # Update XYZu and Bias
                XYZ_r = XYZ_r + dX[0:3]
                Bias = Bias + dX[3]

            LLA_r = ecef2lla(XYZ_r)  # User position in the geodetic coordinate

            # DOP computation
            Lat = (pi * LLA_r[0]) / 180
            Lon = (pi * LLA_r[2]) / 180
            Q = linalg.inv(dot(A.T, A))
            RR = array([[-sin(Lon.real), cos(Lon), [0]],
                        [-sin(Lat) * cos(Lon), -sin(Lat) * sin(Lon), cos(Lat)],
                        [cos(Lat) * cos(Lon), cos(Lat) * sin(Lon), sin(Lat)]]).reshape(3, 3)
            CC = dot(dot(linalg.inv(RR), Q[0:3, 0:3]), RR)
            VDOP = sqrt(CC[2, 2])
            HDOP = sqrt(CC[0, 0] + CC[1, 1])

    return XYZ_r.reshape(1, -1), LLA_r.reshape(1, -1), Bias, HDOP, VDOP, Flag


def EstimateIPP(Lat, Lon, Ele, Azi):
    """"
            Calculate Ionospheric Pierce Point (IPP)
             Inputs:
                    Lat         = Latitude (Degree) of receiver
                    Lon         = Longitude (Degree) of receiver
                    Ele         = Elevation angle (Degree) between satellite and receiver
                    Azi         = Azimuth angle (Degree) between satellite and receiver
             Outputs:
                   Lat_ipp      = Latitude (Degree) of IPP
                   Lon_ipp      = Longitude (Degree) of IPP
                   Heig_ipp     = Height (m) of IPP
            LAB of ThaiPBS, Thailand.
            Version 1 by Somkit Sophan (June 2025)
        """
    Re = 6371230        # The mean radius of the Earth
    h = 350 * (10**3)   # Ionospheric height for mapping function (350 km)
    Z = 90 - Ele
    Z_parm = rad2deg(arcsin((Re * cos(deg2rad(Ele)))/(Re + h)))
    Alpha = Z - Z_parm

    Lat_ipp = rad2deg(arcsin(cos(deg2rad(Alpha)) * sin(deg2rad(Lat)) + sin(deg2rad(Alpha)) * cos(deg2rad(Lat)) * cos(deg2rad(Azi))))

    if ((Lat>70) & (tan(deg2rad(Alpha))*cos(deg2rad(Azi))>tan(deg2rad(90-Lat)))) |\
            ((Lat<-70) & (arctan(deg2rad(Alpha))*cos(deg2rad(Azi))>tan(deg2rad(90-Lat)))):
        Lon_ipp = Lon + 180 - rad2deg(arcsin(sin(deg2rad(Alpha)) * sin(deg2rad(Azi)) / cos(deg2rad(Lat_ipp))))
    else:
        Lon_ipp = Lon + rad2deg(arcsin(sin(deg2rad(Alpha)) * sin(deg2rad(Azi)) / cos(deg2rad(Lat_ipp))))

    Heig_ipp = h

    return [Lat_ipp, Lon_ipp, Heig_ipp]