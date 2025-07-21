from numpy import array, arctan2, sqrt, sin, cos, linalg, pi, radians, dot


def ecef2lla(xyz):
    """
        Convert earth-centered earth-fixed (ECEF)
        cartesian coordinates to latitude, longitude, and altitude
         Input:
            xyz = ECEF coordinates size 1x3 (meter)
         Output:
            lla = Geodetic coordinates size 1x3 (degree)
        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021)
    """

    xyz = xyz.reshape(-1, 1)
    x = xyz[0]
    y = xyz[1]
    z = xyz[2]
    # WGS84 ellipsoid constants:
    a = 6378137  # major-radius of the earth reference ellipsoid
    f = 1 / 298.257223563  # flatting of the earth reference ellipsoid spheroid properties
    b = (1 - f) * a  # Semi-minor axis
    e2 = f * (2 - f)  # Square of (first) eccentricity
    ep2 = e2 / (1 - e2)  # Square of second eccentricity
    pr = sqrt((x ** 2) + (y ** 2))
    beta = arctan2(z, (1 - f) * pr)
    lat = arctan2(z + b * e2 * sin(beta) ** 3, pr - a * e2 * cos(beta) ** 3)
    beta_new = arctan2((1 - f) * sin(lat), cos(lat))
    count = 0
    while count < 5:
        beta = beta_new
        lat = arctan2(z + b * ep2 * sin(beta) ** 3, pr - a * e2 * cos(beta) ** 3)
        beta_new = arctan2((1 - f) * sin(lat), cos(lat))
        count = count + 1

    n = a / sqrt((1 - e2 * sin(lat) ** 2))

    return array([lat * 180 / pi, arctan2(y, x) * 180 / pi, pr * cos(lat) + (z + e2 * n * sin(lat)) * sin(lat) - n])


def lla2ecef(lla):
    """
        Convert geodetic coordinates
        latitude, longitude, and altitude to cartesian coordinates
         Input:
            lla = Geodetic coordinates size 1x3 (degree)
         Output:
            xyz = ECEF coordinates size 1x3 (meter)
        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021)
    """

    # WGS84 ellipsoid constants:
    a = 6378137  # major-radius of the earth reference ellipsoid
    f = 1 / 298.257223563  # flatting of the earth reference ellipsoid Spheroid properties
    e2 = f * (2 - f)  # Square of (first) eccentricity
    lla = lla.reshape(-1, 1)
    lat = lla[0]
    lon = lla[1]
    h = lla[2]
    v = a / sqrt(1 - e2 * sin(radians(lat)) ** 2)
    x = (v + h) * cos(radians(lat)) * cos(radians(lon))
    y = (v + h) * cos(radians(lat)) * sin(radians(lon))
    z = (v * (1 - e2) + h) * sin(radians(lat))
    return array([x, y, z])


def ecef2enu(XYZ_r, XYZ_s):
    """
        Convert earth-centered earth-fixed (ECEF)
        cartesian coordinates to East, North, Up
         Input:
            XYZ_r = Receiver positions size 1x3 (meter)
            XYZ_s = Satellite positions size 1x3 (meter)
         Output:
            enu = Local geodetic coordinates size 1x3 (meter)
        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021)
    """

    XYZ_r = XYZ_r.reshape(-1, 1)
    XYZ_s = XYZ_s.reshape(-1, 1)
    lat, lon, h = ecef2lla(XYZ_r)
    # --- The rotation matrix R of the ECEF coordinates to the local coordinates ---
    e = array([-sin(radians(lon)), cos(radians(lon)), array([0])])
    n = array([-sin(radians(lat)) * cos(radians(lon)), - sin(radians(lat)) * sin(radians(lon)), cos(radians(lat))])
    u = array([cos(radians(lat)) * cos(radians(lon)), cos(radians(lat)) * sin(radians(lon)), sin(radians(lat))])
    r = XYZ_s - XYZ_r
    a = array(r / linalg.norm(r)).reshape(1, -1)  # Unit vector
    # --- The ENU coordinates ---
    east = dot(a, e)
    north = dot(a, n)
    up = dot(a, u)
    return array([east, north, up])


def Rotation_Z(XYZ_s, dT):
    """
        Rotate X & Y axis around the Z axis
         Input:
            XYZ_s = Satellite position in ECEF coordinates size 1x3 (meter)
            dT    = Transmitted time of satellites (sec)
         Output:
            New_XYZ_s = geodetic coordinates size 1x3 (degree)
        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021)
    """
    We = 7.2921151467 * 10 ** -5  # Earth rotation rate(rad/sec)
    Theta = We * dT
    Rz = array([[cos(Theta), sin(Theta), array([0])],
                [-sin(Theta), cos(Theta), array([0])],
                [array([0]), array([0]), array([1])]])
    return dot(XYZ_s, Rz).reshape(1, -1)


def lldistkm(latlon1, latlon2):
    """
        d1km: distance in km based on Haversine formula
        (Haversine: http://en.wikipedia.org/wiki/Haversine_formula)
        d2km: distance in km based on Pythagoras theorem
        (see: http://en.wikipedia.org/wiki/Pythagorean_theorem)
        After: http://www.movable-type.co.uk/scripts/latlong.html

        Inputs:
            latlon1: latlon of origin point [lat lon]
            latlon2: latlon of destination point [lat lon]

        Outputs:
            d1km: distance calculated by Haversine formula
            d2km: distance calculated based on Pythagoras theorem

        Example 1, short distance:
            latlon1=[-43 172];
            latlon2=[-44  171];
            [d1km d2km]=distance(latlon1,latlon2)
            d1km =  137.365669065197 (km)
            d2km =  137.368179013869 (km)

            d1km approximately equal to d2km

        Example 2, longer distance:
             latlon1=[-43 172];
             latlon2=[20  -108];
             [d1km d2km]=distance(latlon1,latlon2)
             d1km = 10734.8931427602 (km)
             d2km = 31303.4535270825 (km)
         d1km is significantly different from d2km (d2km is not able to work for longer distances).

        CSSRG-LAB of KMITL, Thailand.
        Version 1 by Somkit Sophan (July 2021)
    """
    radius = 6371
    lat1 = latlon1[0] * pi / 180
    lat2 = latlon2[0]*pi/180
    lon1 = latlon1[1]*pi/180
    lon2 = latlon2[1]*pi/180
    deltaLat = lat2-lat1
    deltaLon = lon2-lon1
    a = sin(deltaLat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(deltaLon / 2) ** 2
    c = 2 * arctan2(sqrt(a), sqrt(1 - a))
    d1km = radius*c     # Haversine distance

    # x = deltaLon * cos((lat1 + lat2) / 2)
    # y = deltaLat
    # d2km = radius * sqrt(x * x + y * y)   # Pythagoras distance

    return d1km
