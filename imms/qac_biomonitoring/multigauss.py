import numpy as np
from numpy import argmax
import pandas as pd
import os
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, gaussian, convolve
from scipy.optimize import curve_fit

# Smooth the raw EIC with a Gaussian convolution
def process_data(rt, rt_i, file_name, mz, sample_type): # Raw values extracted usisng dhrmasslynxapi 
    window_len = 51
    window = gaussian(window_len, std=1) # Adjust to change degree of smoothing 
    smoothed_intensity = convolve(rt_i, window/window.sum(), mode="same")

    # Identify peaks in the smoothed data
    peak_indices, _ = find_peaks(smoothed_intensity, prominence=1500, distance=5) # Adjust peak picking filters

    # Add identified peaks to rt_list for output file
    rt_list = []
    for j in peak_indices:
        label = float(rt[j])
        round_label = str(round(label, 2))
        rt_list.append(round_label)

    # Extract mu (mean) values for the fixed Gaussian functions from the smoothed data
    mu_values = rt[peak_indices]

    # Initial guesses for amplitude (A) and sigma from the smoothed data
    A_guesses = smoothed_intensity[peak_indices]
    sigma_guesses = [0.05] * len(peak_indices)
    p0 = [val for sublist in zip(A_guesses, sigma_guesses) for val in sublist]

    # Define the multi-Gaussian function with fixed mu values
    def multi_gaussian_fixed_mu(x, *params):
        y = np.zeros_like(x)
        for i in range(0, len(params), 2):
            A = params[i]
            sigma = params[i + 1]
            y += A * np.exp(-((x - mu_values[int(i/2)])**2) / (2 * sigma**2))
        return y

    # Create "Extracted Chromatograms" folder if it does not exist
    chromatogram_directory = "Extracted Chromatograms"
    if not os.path.exists(chromatogram_directory):
        os.makedirs(chromatogram_directory)

    # Attempt to fit the smoothed data to a multi-Gaussian function           
    try:
        popt_multi_gaussian, _ = curve_fit(multi_gaussian_fixed_mu, rt, rt_i, p0=p0, maxfev=20000)

        areas = []
        for i in range(0, len(popt_multi_gaussian), 2):
            A = popt_multi_gaussian[i]
            sigma = popt_multi_gaussian[i + 1]
            mu = mu_values[int(i/2)]
            area = A * sigma * np.sqrt(2 * np.pi)  # Gaussian integral
            areas.append(area)

        # Set up LC chromatogram figure for successfully fitted EICs
        title_lc = "Extracted Chromatogram ({}) \nm/z: {:.4f} \u00B1 0.025".format(file_name, mz)
        fname_lc = "{}/{}_{}_{}_EIC.png".format(chromatogram_directory, file_name, mz, sample_type.capitalize())
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        ax.plot(rt, rt_i, "b-", lw=1.5, label="Raw Data") # Plot raw EIC
        ax.plot(rt, smoothed_intensity, "g--", lw=1.5, label="Smoothed Data") # Plot smoothed data
        ax.plot(rt, multi_gaussian_fixed_mu(rt, *popt_multi_gaussian), "r-", lw=1.5, label="Gaussian Fit") # Plot fitted data
        y_values = multi_gaussian_fixed_mu(rt[peak_indices], *popt_multi_gaussian)
        for j in peak_indices:
            label = str(rt[j])
            plt.annotate(label[:4], xy=(rt[j], rt_i[j]), xytext=(rt[j]+0.04, rt_i[j] * 0.95), fontsize=10, fontweight="bold", fontname="Arial") # Annotate rt
        plt.scatter(rt[peak_indices], y_values, color="purple", marker="*", s=40, label="Identified Peaks") # Mark identified peaks with a star
        title_fontprops = {"fontsize": 12, "fontweight": "bold", "fontname": "Arial"} 
        axes_label_fontprops = {"fontsize": 12, "fontweight": "bold", "fontname": "Arial"}
        tick_label_fontprops = {"weight": "bold", "family": "Arial", "size": 10}
        ax.set_title(title_lc, **title_fontprops)
        ax.set_xlabel("Retention Time [min]", **axes_label_fontprops)
        ax.set_ylabel("Intensity", **axes_label_fontprops)  
        legend = ax.legend(loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white") # Create legend
        for text in legend.get_texts():
            text.set_fontname("Arial")
        ax.tick_params(axis="both", which="major", labelsize=10, width=1)
        for label in ax.get_xticklabels():
            label.set_fontname("Arial")
            label.set_fontweight("bold")
            label.set_fontsize(10)
        max_intensity_y_limit = max(rt_i) + 0.1 * max(rt_i)
        y_tick_values = np.linspace(0, max_intensity_y_limit, num=10, endpoint=True)  # Generate 10 evenly spaced tick values
        ax.set_yticks(y_tick_values)
        ax.set_yticklabels([int(y) for y in y_tick_values], fontdict=tick_label_fontprops)
        ax.set_ylim(0, max_intensity_y_limit) # Modify y axis range
        plt.xlim(rt[peak_indices[0]] - 1, rt[peak_indices[-1]] + 1) # Modify x axis range
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1.5)
        ax.spines["bottom"].set_linewidth(1.5)
        plt.tight_layout()
        plt.savefig(fname_lc, dpi=350, bbox_inches="tight")
        plt.close()

    # Set up LC chromatogram figure for data that cannot be fit 
    except RuntimeError:
        title_lc = "Extracted Chromatogram ({}) \nm/z: {:.4f} \u00B1 0.025".format(file_name, mz)
        fname_lc = "{}/{}_{}_{}_EIC.png".format(chromatogram_directory, file_name, mz, sample_type.capitalize())
        fig, ax = plt.subplots(figsize=(6.4, 4.8))
        fig = plt.figure() # Generate LC chromatogram figure with only raw data
        ax = fig.add_subplot(111)
        ax.plot(rt, rt_i, "b-", lw=1.5, label="Raw Data")
        title_fontprops = {"fontsize": 12, "fontweight": "bold", "fontname": "Arial"} 
        axes_label_fontprops = {"fontsize": 12, "fontweight": "bold", "fontname": "Arial"}
        tick_label_fontprops = {"weight": "bold", "family": "Arial", "size": 10}  
        ax.set_title(title_lc, **title_fontprops)
        ax.set_xlabel("Retention Time [min]", **axes_label_fontprops)
        ax.set_ylabel("Intensity", **axes_label_fontprops)  
        legend = ax.legend(loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white")
        for text in legend.get_texts():
            text.set_fontname("Arial")
        legend.get_frame().set_edgecolor("black")
        legend.get_frame().set_linewidth(1)
        ax.tick_params(axis="both", which="major", labelsize=8, width=0.5)
        ax.set_xticks(ax.get_xticks()) # Format tick labels manually
        ax.set_xticklabels(ax.get_xticks(), fontdict=tick_label_fontprops)
        max_intensity = max(rt_i)
        y_tick_values = np.linspace(0, max_intensity + 1000, num=10, endpoint=True)  # Generate 10 evenly spaced tick values
        ax.set_yticks(y_tick_values)
        ax.set_yticklabels([int(y) for y in y_tick_values], fontdict=tick_label_fontprops)
        plt.ylim(0, max_intensity + 1000) # Modify y axis range
        max_rt = rt[argmax(rt_i)]
        plt.xlim(max_rt - 4, max_rt + 4) # Modify x axis range
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1.5)
        ax.spines["bottom"].set_linewidth(1.5)
        plt.tight_layout()
        plt.savefig(fname_lc, dpi=350, bbox_inches="tight")
        plt.close() 
        return [], [], []

    return peak_indices, rt_list, areas

