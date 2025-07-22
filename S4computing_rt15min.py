import os
import glob
import datetime as dt
import time
from matplotlib import pyplot   # Install packages
from numpy import array, ones, mean, nanmean, append, abs, fix, real, sqrt, isnan, nan, delete  # Install packages
from pandas import DataFrame    # Install packages
from gnsspy.io.readFile import read_obsFile, read_navFile
from Functions.Coordinate import ecef2lla
from Functions.EstimatePosition import GPS_xyz, EstimateIPP

# --- Define PRN of GPS satellites ---
PRN_list = ['G01', 'G02', 'G03', 'G04', 'G05', 'G06', 'G07', 'G08', 'G09', 'G10',
            'G11', 'G12', 'G13', 'G14', 'G15', 'G16', 'G17', 'G18', 'G19', 'G20',
            'G21', 'G22', 'G23', 'G24', 'G25', 'G26', 'G27', 'G28', 'G29', 'G30',
            'G31', 'G32']

# === Initial setup ===
Input_path = "C:/S4C_Computations/Results/"  # *** Input path ***
Output_path = "C:/S4C_Computations/Results/"   # *** Output path ***
CNR_mask = 30    # C/N0 mask
Ele_mask = 15    # Elevation mask
Station = "SNxxx"   # Station code
Satellite = "GPS"   # Satellite system
Windows_SI_det = 10     # Window size for computing the SI_det values
Windows_S4C = 10    # Window size for computing the S4C values
Windows_Sm = 15     # Window size for smoothing the S4C values

# === Begin positioning error calculation ===
for root, dirs, files in os.walk(Input_path):
    if len(root):
        Obs_files = glob.glob(root + Station + "*.*o", recursive=True)
        Nav_files = glob.glob(root + Station + "*.*n", recursive=True)
        if len(Obs_files) & len(Nav_files):
            for i in range(0, len(Obs_files), 1):
                if Obs_files[i][-8:-1] == Nav_files[i][-8:-1]:
                    Obs_name = Obs_files[i]  # *** Navigation file name
                    Nav_name = Nav_files[i]  # *** Observation file name
                else:
                    continue

                # --- Read RINEX files ---
                Nav = read_navFile(Nav_name)  # Read navigation file
                Obs = read_obsFile(Obs_name)  # Read observation file

                XYZ_Obs = array(Obs.approx_position)  # Reference position (ECEF)
                LLA_Obs = ecef2lla(XYZ_Obs)   # Convert the ECEF to [Lat, Lon, Height]

                # --- Set sampling time ---
                Epoch = dt.datetime(Obs.time_first[0], Obs.time_first[1], Obs.time_first[2],
                                    Obs.time_first[3], Obs.time_first[4], Obs.time_first[5])

                Obs_SYS = Obs.observation['SYSTEM'] == Satellite # Select GPS
                Obs_C1 = Obs.observation.C1C[Obs_SYS]  # Read code-pseudorange L1
                Obs_S1 = Obs.observation.S1C[Obs_SYS]  # Read C/N0 L1

                Epoch_temp = Obs_C1.index.levels[0]

                S1 = DataFrame(ones((1, PRN_list.__len__())) * nan,
                               index=[Epoch_temp[0]], columns=PRN_list)     # L1 C/N0
                SI_det = DataFrame(ones((1, PRN_list.__len__())) * nan,
                                   index=[Epoch_temp[0]+dt.timedelta(seconds=Windows_SI_det-1)], columns=PRN_list) # SI_det of L1 C/N0
                S4C_L1 = DataFrame(ones((1, PRN_list.__len__())) * nan,
                                   index=[Epoch_temp[0]+dt.timedelta(seconds=Windows_SI_det+Windows_S4C-1)], columns=PRN_list)
                S4C_L1_Sm = DataFrame(ones((1, PRN_list.__len__())) * nan,
                                   index=[Epoch_temp[0] + dt.timedelta(seconds=Windows_SI_det+Windows_S4C-1)], columns=PRN_list)
                Lat_ipp = DataFrame(ones((1, PRN_list.__len__())) * nan,
                                    index=[Epoch_temp[0] + dt.timedelta(seconds=Windows_SI_det + Windows_S4C - 1)], columns=PRN_list)
                Lon_ipp = DataFrame(ones((1, PRN_list.__len__())) * nan,
                                    index=[Epoch_temp[0] + dt.timedelta(seconds=Windows_SI_det + Windows_S4C - 1)], columns=PRN_list)

                Start = time.time()
                Count = 1
                print("Compute S4Cs...")
                for ind_t in Obs_C1.index.levels[0]:  # A sampling time loop
                    try:
                        Ci = Obs_C1[ind_t]
                        Si = Obs_S1[ind_t]
                        PRN = Obs_C1[ind_t].index.values
                        SOD = ind_t.hour * 3600 + ind_t.minute * 60 + ind_t.second  # Second of day
                    except:
                        continue

                    ind_n = isnan(Ci.values) | isnan(Si.values)  # Matched information between Ci and Si
                    PRN = delete(PRN, ind_n)

                    S1.loc[ind_t, :] = nan
                    if S1.__len__() >= Windows_SI_det:
                        SI_det.loc[ind_t, :] = nan
                    if (SI_det.__len__() > Windows_S4C) & ((Count % Windows_S4C) == 0):
                        S4C_L1.loc[ind_t, :] = nan
                        if S4C_L1.__len__() >= Windows_Sm:
                            S4C_L1_Sm.loc[ind_t, :] = nan
                            Lat_ipp.loc[ind_t, :] = nan
                            Lon_ipp.loc[ind_t, :] = nan

                    for ind_i in PRN:  # A PRN loop
                        # Satellite position estimations
                        if not isnan(Ci[ind_i]):
                            LLA_ipp = [nan, nan, nan]
                            try:
                                XYZ_s, Clock_s, Ele_temp, Azi_temp = \
                                    GPS_xyz(Nav.ephemeris, Epoch, SOD, Ci[ind_i], ind_i, XYZ_Obs)    # Estimated function
                                if (Ele_temp[0][0] >= Ele_mask) & (Si[ind_i] >= CNR_mask):
                                    S1.loc[ind_t,ind_i] = Si[ind_i]           # C/N0 log
                                    if S4C_L1.__len__() >= Windows_Sm:
                                        LLA_ipp = EstimateIPP(LLA_Obs[0][0], LLA_Obs[1][0], Ele_temp[0][0], Azi_temp[0][0])
                            except:
                                continue

                            # Compute SI_det
                            if S1.__len__() >= Windows_SI_det:
                                S1_N0 = 10 ** (0.1 * S1)  # Convert dB to Watt
                                SI_det.loc[ind_t, ind_i] = S1_N0[ind_i][ind_t] / sum(S1_N0.loc[:, ind_i])

                            # Compute S4C values
                            if (SI_det.__len__() > Windows_S4C) & ((Count % Windows_S4C) == 0):
                                Temp = SI_det.loc[:(ind_t-dt.timedelta(seconds=1)), ind_i]
                                if sum(isnan(Temp)) == 0:
                                    S4C_L1.loc[ind_t, ind_i] = real(sqrt((mean(Temp ** 2) - (mean(Temp) ** 2)) / (mean(Temp) ** 2))) / 0.5

                                    # Smooth S4C values
                                    if S4C_L1.__len__() >= Windows_Sm:
                                        S4C_L1_Sm.loc[ind_t, ind_i] = mean(S4C_L1.loc[:, ind_i])
                                        Lat_ipp.loc[ind_t, ind_i] = LLA_ipp[0]
                                        Lon_ipp.loc[ind_t, ind_i] = LLA_ipp[1]
                        else:
                            continue

                    if S1.__len__() >= Windows_SI_det:
                        S1 = S1.drop(S1.index[0])

                    if SI_det.__len__() > Windows_S4C:
                        SI_det = SI_det.drop(SI_det.index[0])

                    if S4C_L1.__len__() >= Windows_Sm:
                        S4C_L1 = S4C_L1.drop(S4C_L1.index[0])

                    Count += 1

                Finish = time.time()
                print("Compute S4Cs in", "{0:.2f}".format(Finish - Start), "seconds.")

                # --- Save S4C calculation ---
                S4C_L1_Sm.to_csv(Output_path + Station + "_S4C_last15min.csv")
                Lat_ipp.to_csv(Output_path + Station + "_Lat_last15min.csv")
                Lon_ipp.to_csv(Output_path + Station + "_Lon_last15min.csv")

                # # --- plot ---
                # pyplot.cla()
                # pyplot.plot(S4C_L1_Sm)
                # pyplot.xlabel('UTC')
                # pyplot.ylabel('S4C')
                # pyplot.ylim((0, 1))
                # pyplot.grid(True)
                # pyplot.show()