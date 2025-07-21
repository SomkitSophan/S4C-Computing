from numpy import isnan, sum, delete, nan
import datetime as dt


def CarrierSmoothedCode(epoch, code, carrier, pre_code, pre_carrier, lambda_carrier, sec):
    """
    Objective: To smooth the code pseudorange by using the carrier pseudorange.
    Input:
            epoch       = epoch time data
            code        = code pseudorange data
            carrier     = carrier-phase pseudorange data
            pre_code    = code (t-1)
            pre_carrier = carrier (t-1)
            lambda_carrier  = Wavelength of carrier-phase
            sec         = sampling time
    Output:
            code_smooth = carrier smoothed code data
            PRN         = PRNs of the GNSS satellite
            pre_code    = code pseudorange data for (t+1)
            pre_carrier = carrier-phase pseudorange data for (t+1)
    CSSRG-LAB of KMITL, Thailand.
    Version 1 by Somkit Sophan (July 2021)
    """
    # Define smoothed values
    a = 100

    # Find data by the sampling time "sec"
    ind_t = epoch + dt.timedelta(seconds=sec)
    if sum(code.index.levels[0] == ind_t) > 0:
        Ci = code[ind_t]
        Li = carrier[ind_t]
        PRN = code[ind_t].index.values
        # print(sec)  # Print the sampling time number
    else:
        Flag = 0
        code_smooth = nan
        Ci = nan
        PRN = nan
        # print("Data is empty")
        return code_smooth, Ci, PRN, pre_code, pre_carrier, Flag

    Ind_n = isnan(Ci.values) | isnan(Li.values)     # Matched information between Ci and Li
    PRN = delete(PRN, Ind_n)

    # === The carrier phase smoothed method ===
    if sum(pre_code[0][PRN].values == 0) > 0:   #
        Ind = pre_code[0][PRN].values == 0
        pre_code.loc[PRN[Ind],0] = Ci[PRN[Ind]].values
        pre_carrier.loc[PRN[Ind],0] = Li[PRN[Ind]].values

    if len(PRN) > 3:  # Need at least 4 satellites to estimate the receiver positions
        code_smooth = (1 / a) * Ci[PRN].values + ((a - 1) / a) * (
                pre_code[0][PRN].values + lambda_carrier * (Li[PRN].values - pre_carrier[0][PRN].values))
        pre_code.loc[PRN,0] = code_smooth
        pre_carrier.loc[PRN,0] = Li[PRN].values
        Flag = 1
        return code_smooth.reshape(-1, 1), Ci[PRN].values.reshape(-1, 1), PRN, pre_code, pre_carrier, Flag
    else:
        Flag = 0
        code_smooth = nan
        # print("GPS data less than 4 satellites")
        return code_smooth, Ci, PRN, pre_code, pre_carrier, Flag
