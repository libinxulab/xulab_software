import pandas as pd
import numpy as np
from numpy.polynomial.polynomial import Polynomial
from numpy import argmax, sqrt
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, gaussian, convolve
from scipy import asarray as ar, exp
import math
from matplotlib import pyplot as plt
from matplotlib import rcParams
import os
from dhrmasslynxapi.reader import MassLynxReader
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRawXL
from multigauss import process_data
from filter import filter_data

# Basic Gaussian function for peak fitting
def gauss(x, A, B, C):
    if abs(C) < 0.01:
        C = 0.01
    return A * np.exp(-(x - B) ** 2 / (2 * C ** 2))

# Fit the extracted mobilogram with a Gaussian function and return the fitted parameters
def peak_fit(dt, dt_i, p0="guess"):  # p0 are initial parameters for peak fit
    if p0 == "guess":
        p0 = (max(dt_i), dt[argmax(dt_i)], 0.5)
    opt, cov = curve_fit(gauss, dt, dt_i, maxfev=5000, p0=p0)
    return opt

# Set global font conditions for figures
params = {"font.size": 8,
         "font.family": "Arial",
          "font.weight": "bold"}
plt.rcParams.update(params)

# Set up mobilogram figure
def atd(t, t_refined, dt_i, fit_i, A, B, title_atd, fname_atd):
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    ax.plot(t, dt_i, "o--b", lw=1.5, ms=2, label="Raw Data") # Plot raw mobilogram
    ax.plot(t_refined, fit_i, "r-", lw=1.5, label="Gaussian Fit") # Plot fitted data
    legend = ax.legend(loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white") # Create legend 
    for text in legend.get_texts():
        text.set_fontname("Arial")
    ax.text(B + 0.1, 0.95 * max(fit_i), "{:.2f}".format(B), c="k", fontsize=10, fontweight="bold", fontname="Arial") # Annotate peak of fitted data (dt)
    ax.set_title(title_atd, fontsize=12, fontweight="bold", fontname="Arial")
    ax.set_xlabel("Drift Time [ms]", fontsize=12, fontweight="bold", fontname="Arial")
    ax.set_ylabel("Intensity", fontsize=12, fontweight="bold", fontname="Arial")
    max_dt_i_y_limit = max(dt_i) + 0.1 * max(dt_i)
    plt.ylim(0, max_dt_i_y_limit) # Modify y axis range
    plt.xlim(0, B + 2) # Modify x axis range
    tick_label_fontprops = {"weight": "bold", "family": "Arial", "fontsize": 10}
    ax.tick_params(axis="both", which="major", labelsize=10, width=1)
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xticklabels(ax.get_xticks(), fontdict=tick_label_fontprops)
    y_ticks = ax.get_yticks()
    ax.set_yticklabels([int(y) for y in y_ticks], fontdict=tick_label_fontprops)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    plt.tight_layout()
    plt.savefig(fname_atd, dpi=350, bbox_inches="tight")
    plt.close()

output_rows = []

# Main execution sequence
def main(input_file, calibration_file):
    df_input = pd.read_excel(input_file) # Read input file into a pandas dataframe
    df_input.reset_index(drop=True, inplace=True)

    for i, row in df_input.iterrows():

        # Extract the values for file_name, m/z, and sample_type
        file_name = row["file_name"]
        mz = row["mz"]
        sample_type = row["sample_type"]

        # Print the current m/z being queried to the terminal
        print("\n\traw file: {} m/z: {}".format(file_name, mz))

        # Extract mobilogram from the raw file using m/z target
        rdr = MassLynxReader(file_name)
        mz = float(mz)
        t, dt_i = rdr.get_chrom(2, float(mz), 0.025)  # Adjust the MS function and mass tolerance window

        # Parameter for smoothing Gaussian curve
        t_refined = np.arange(min(t), max(t), float(0.01))

        # Initialize Gaussian function
        A, B, C = peak_fit(t, dt_i)

        # Fit the raw ATD to Gaussian function
        fit_i = gauss(t_refined, A, B, C)

        # Round fitted parameter B (which equals dt)
        round_B = str(round(B, 2))

        # Apply FWHM and intensity thresholds for extracted mobilogram
        fwhm = C * 2.355
        if not (fwhm < 2.5 and fwhm > 0.05 and max(dt_i) > 500): # Adjust FWHM and ATD intensity filters
            round_B = None

        # Convert dt to float
        dt = float(round_B) if round_B else None

        # Calculate calibrated CCS value from calibration curve
        ccs = None
        if pd.notnull(dt):
            ccs = cal_data.calibrated_ccs(mz, dt)
            ccs = round(ccs, 2)

        # Extract m/z-selected EIC from the raw data
        rt, rt_i = rdr.get_chrom(0, float(mz), 0.025)  # Adjust the MS function and mass tolerance window

        # Normalize the EIC for optional filtering and/or generating figures; not currently utilized
        max_intensity = max(rt_i)
        normalized_i = [intensity / max_intensity * 100 for intensity in rt_i]

        # Store rt and rt_i as numpy arrays for processing
        rt = np.array(rt)
        rt_i = np.array(rt_i)

        # Smooth and fit data with multi-Gaussian function 
        data = (rt, rt_i, file_name, mz, sample_type)
        peak_indices, rt_list, areas = process_data(*data)

        # Append the rt and peak areas to the output_rows list
        for rt_value, peak_area in zip(rt_list, areas):
            output_rows.append([file_name, mz, sample_type, row["ccs_calibrant"], row["gradient"], row["column_type"], dt, ccs, rt_value, peak_area])

        # Create title for mobilogram figure
        title_atd = "Extracted Mobilogram ({}) \nm/z: {:.4f} \u00B1 0.025 → FWHM ~ {:.2f} ms".format(file_name, mz, fwhm)

        # Create file name for mobilogram figure
        fname_atd = "{}/{}_{}_{}_EIM.png".format(mobilogram_directory, file_name, mz, sample_type.capitalize()) # Change where mobilogram figures are saved

        # Generate mobilogram figure 
        atd(t, t_refined, dt_i, fit_i, A, B, title_atd, fname_atd) # Comment out this code if figures are not needed

    # Create a new dataframe from the output_rows list
    df_output = pd.DataFrame(output_rows, columns=["file_name", "mz", "sample_type", "ccs_calibrant", "gradient", "column_type", "dt", "ccs", "rt", "peak_area"])

    # Write the dataframe to an Excel file
    output_file = input_file[:-5] + "_processed.xlsx"
    df_output.to_excel(output_file, index=False)
    print("\nProcessing successful. View results in {}".format(output_file))

    # Filter sample data based on corresponding control peak areas
    filtered_output_file = filter_data(output_file)
    print("\nFiltering successful. View results in {}\n".format(filtered_output_file))

if __name__ == "__main__":
    input_file = "2023_07_28_RO1_feces_30B_100B.xlsx" # Enter path to the processing template containing file names
    calibration_file = "2023_07_28_polyala_mz_ccs.xlsx" # Enter path to the CCS calibrant 
    cal_data = CCSCalibrationRawXL(calibration_file) # Initialize CCS calibration object and load calibrant Excel spreadsheet
    curve_fname = calibration_file[:-5]
    cal_data.cal_curve_figure(curve_fname) # Generate a .png image of the CCS calibration curve
    print("\n...Calibration successful. Extracting mobilograms, calibrated CCS values, and LC chromatograms from raw data...")
    mobilogram_directory = "Extracted Mobilograms" # Create "Extracted Mobilograms" if it does not exist
    if not os.path.exists(mobilogram_directory):
        os.makedirs(mobilogram_directory)
    main(input_file, calibration_file)