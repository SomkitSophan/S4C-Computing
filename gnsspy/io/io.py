"""
Class definitions for I/O opreations
"""
# ======================================================================
class Observation:
    """
    Observations class for RINEX Observation (*.*o) files
    """
    def __init__(self, filename=None, epoch=None, observation=None, approx_position=None,
                 receiver_type=None, antenna_type=None, interval=None,
                 receiver_clock=None, version=None, observation_types=None,
                 time_first=None, time_last=None):
        self.filename          = filename 
        self.epoch             = epoch 
        self.observation       = observation
        self.approx_position   = approx_position
        self.receiver_type     = receiver_type
        self.antenna_type      = antenna_type
        self.interval          = interval
        self.receiver_clock    = receiver_clock
        self.version           = version
        self.observation_types = observation_types
        self.time_first        = time_first
        self.time_last         = time_last
# ======================================================================
class _ObservationTypes:
    def __init__(self, ToB_GPS=None, ToB_GLONASS=None, ToB_GALILEO=None,
                 ToB_COMPASS=None, ToB_QZSS=None, ToB_IRSS=None, ToB_SBAS=None):
        self.GPS     = ToB_GPS
        self.GLONASS = ToB_GLONASS
        self.GALILEO = ToB_GALILEO
        self.COMPASS = ToB_COMPASS
        self.QZSS    = ToB_QZSS
        self.IRSS    = ToB_IRSS
        self.SBAS    = ToB_SBAS
# ======================================================================
class Navigation:
    """
    Navigation class for RINEX Observation (*.*n/p) files
    """
    def __init__(self, ephemeris=None, Ion=None, dUTC=None, version=None):
        self.ephemeris = ephemeris
        self.Ion = Ion
        self.dUTC = dUTC
        self.version = version
# ======================================================================
class Navigation_DEPRECATED:
    """
    Broadcast Ephemeris in RINEX file
    """
    def __init__(self, PRN=None, year=None, month=None, day=None, hour=None,
                 min=None, sec=None, af0=None, af1=None, af2=None, IODE=None,
                  crs=None, delta_n=None, M0=None, cuc=None, e=None, cus=None,
                 sqrtA=None, Toe=None, cic=None, OMEGA0=None, cis=None, i0=None,
                 crc=None, omega=None, OMEGAdot=None, Idot=None, code2=None,
                 SOW=None, flagP2=None, SVa=None, SVh=None, TGD=None, IODC=None,
                 Ts=None, dT=None):
        # ---
        self.PRN = PRN
        self.year = year
        self.month = month
        self.day = day
        self.hour = hour
        self.min = min
        self.sec = sec
        self.af0 = af0
        self.af1 = af1
        self.af2 = af2
        # ---
        self.IODE = IODE
        self.crs = crs
        self.delta_n = delta_n
        self.M0 = M0
        # ---
        self.cuc = cuc
        self.e = e
        self.cus = cus
        self.sqrtA = sqrtA
        # ---
        self.Toe = Toe
        self.cic = cic
        self.OMEGA0 = OMEGA0
        self.cis = cis
        # ---
        self.i0 = i0
        self.crc = crc
        self.omega = omega
        self.OMEGAdot = OMEGAdot
        # ---
        self.Idot = Idot
        self.code2 = code2
        self.SOW = SOW
        self.flagP2 = flagP2
        # ---
        self.SVa = SVa
        self.SVh = SVh
        self.TGD = TGD
        self.IODC = IODC
        # ---
        self.Ts = Ts
        self.dT = dT
# ======================================================================
class PEphemeris:
    """
    Class definition for SP3 file (Precise Ephemeris)
    """
    def __init__(self, epoch=None, ephemeris=None):
        self.epoch = epoch
        self.ephemeris = ephemeris
