import numpy as np
import pandas as pd
from inspect import signature


def HST_detrending(lightcurve, v1, v2, a, b, t0, P):
    return


def get_ldtk_params(model):
    if model == "linear":
        return ["q0"]
    if model in ["quadratic", "squareroot", "power2"]:
        return ["q0", "q1"]
    if model == "nonlinear":
        return ["q0", "q1", "q2", "q3"]


def count_number_lcs(
    inputdata,
    other_params=["a", "inc", "rp"],
    limb_darkening_model="quadratic",
    ld_fit_method="coupled",
    normalise=True,
    detrend=True,
    error_scaling=True,
    max_batch_parameters=35,
    detrending_list=[["nth order", 2]],
):
    total_params = 3  # "P", "P_prime", "P_dprime"
    multipliers = {"P": 1, "P_prime": 1, "P_dprime": 1}
    if "a" in other_params:
        total_params += 1
        multipliers["a"] = 1

    if "inc" in other_params:
        total_params += 1
        multipliers["inc"] = 1

    if normalise:
        if not "norm" in other_params:
            other_params += ["norm"]

    if error_scaling:
        if not "escale" in other_params:
            other_params += ["escale"]

    if ld_fit_method == "off":
        ldtk_params_count = 0
        ldtk_params = []

    else:
        ldtk_params = get_ldtk_params(limb_darkening_model)

        ldtk_params_count = len(ldtk_params)

    df = pd.read_csv(inputdata)
    df = df.sort_values("Epochs", ignore_index=True)

    total_params_list = np.empty(0)

    num_detrending_params = np.empty(0)
    if detrend:
        for i, j in enumerate(detrending_list):
            try:
                num_detrending_params = np.append(num_detrending_params, int(j[1]))
            except:
                sig = signature(j[1])
                params = sig.parameters

                # First parameter is lightcurve, last two parameters are t0 and P
                num_detrending_params = np.append(
                    num_detrending_params, len(params) - 3
                )

        detrending_coeffs = []
        for j, val in enumerate(num_detrending_params):
            for i in range(int(val)):
                detrending_coeffs.append("d" + str(j) + "_" + str(i))
    else:
        num_detrending_params = 0
        detrending_coeffs = []

    other_params = other_params + ldtk_params + detrending_coeffs

    for i in range(2, 35):
        _total_params = 0 + total_params

        indices = np.unique(np.linspace(0, len(df) - 1, i, dtype=int))
        selected_df = df.iloc[indices]

        num_lightcurves = len(selected_df)
        num_filters = len(np.unique(selected_df["Filter"].to_numpy(dtype=float)))

        if "rp" in other_params:
            _total_params += num_filters

        _total_params += num_filters * ldtk_params_count

        if normalise:
            _total_params += num_lightcurves

        if error_scaling:
            _total_params += num_lightcurves

        if detrend:
            detrending = selected_df["Detrending"].to_numpy(dtype=float)
            _, counts_detrending = np.unique(detrending, return_counts=True)
            _total_params += sum(num_detrending_params * counts_detrending)

        total_params_list = np.append(total_params_list, int(_total_params))
        if _total_params > max_batch_parameters:
            break

        else:
            if "rp" in other_params:
                multipliers["rp"] = num_filters
            if normalise:
                multipliers["norm"] = num_lightcurves
            if error_scaling:
                multipliers["escale"] = num_lightcurves
            if detrend:
                for j, val in enumerate(num_detrending_params):
                    for i in range(int(val)):
                        multipliers["d" + str(j) + "_" + str(i)] = counts_detrending[j]

            for l in ldtk_params:
                multipliers[l] = num_filters
            number_lcs = num_lightcurves

    print("Parameters fitted: (param:counts)\n", multipliers)
    print("Parameters in this batch:", int(total_params_list[-2]))
    print("Number of lightcurves used:")
    if total_params_list[-1] == total_params_list[0]:
        print(len(df))
        return len(df)

    else:
        print(int(number_lcs))
        return int(number_lcs)
    

     
