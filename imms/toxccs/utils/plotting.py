"""
toxccs/utils/plotting.py

Ryan Nguyen (ryan97@uw.edu)
12/21/24

description:
        Module with utilities for plotting extracting data. Current module treats single and multi-Gaussian extracted ion mobilograms separately.
"""

import numpy as np
import warnings
from matplotlib import pyplot as plt
from scipy.signal import convolve, gaussian
from PIL import Image
from toxccs.utils.smiles import smiles_to_structure
from toxccs.utils.gaussian import (
    process_chromatogram,
    process_multi_chromatogram,
    multi_gaussian_fixed_mu,
)


def combined_figure(
    rt,
    rt_i,
    t,
    dt_i,
    t_refined,
    fit_i,
    A,
    B,
    file_name,
    compound_name,
    adduct,
    mz,
    monoisotopic_mz,
    monoisotopic_intensity,
    rt_value,
    fwhm,
    mu_values,
    fname_combined,
    mz_spectrum,
    intensity_spectrum,
    smiles,
):
    """
    toxccs/utils/plotting.py/combined_figure
    description:
            Generates a figure containing the MS1 scan, extracted ion chromatogram, and extracted ion mobilogram with their Gaussian-fitted curves.
    parameters:
            rt (list or numpy.ndarray) -- array of retention times.
            rt_i (list or numpy.ndarray) -- corresponding intensity values for each retention time.
            t (list or numpy.ndarray) -- raw EIM time points.
            dt_i (list or numpy.ndarray) -- raw EIM intensity values.
            t_refined (list or numpy.ndarray) -- refined data points for plotting the fit.
            fit_i (list or numpy.ndarray) -- intensity values of the fitted curve.
            A (float) -- Gaussian parameter.
            B (float) -- Gaussian parameter.
            file_name (str) --name of the .raw file from which the data is extracted.
            compound_name (str) -- name of the compound being analyzed.
            adduct (str) -- adduct type of the compound being analyzed.
            mz (float) -- target m/z value used to extract the raw data.
            monoisotopic_mz (float) -- observed monoisotopic m/z value.
            rt_value (float) -- observed retention time.
            fwhm (float) -- FWHM of the EIM.
            mu_values -- (ndarray) -- mu values for multi-Gaussian function.
            fname_combined (str) -- file name for saving the plot.
            smoothed_ms1_data (list or numpy.ndarray) -- smoothed MS1 data.
            smiles (str) -- SMILES string of the chemical being processed.
    returns:
            (.png) -- combined image of the extracted and fitted data.
    """

    window_len = 51
    window = gaussian(window_len, std=1)
    smoothed_rt_i = convolve(rt_i, window / window.sum(), mode="same")

    # Retrieve LC chromatogram information
    peak_indices, rt_list, areas, peak_ranges, popt_multi_gaussian, mu_values = (
        process_chromatogram(
            rt,
            rt_i,
            file_name,
            compound_name,
            adduct,
            mz,
            "sample",
            generate_images=False,
        )
    )

    # Set up plot
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(6.4, 9.6))

    # Plot LC chromatogram figure for successfully fitted EICs
    ax2.plot(rt, rt_i, "b-", lw=1.5, label="Raw Data")
    ax2.plot(rt, smoothed_rt_i, "g--", lw=1.5, label="Smoothed Data")

    # Condition to plot if peaks were successfully fitted
    if len(peak_indices) > 0:
        ax2.plot(
            rt,
            multi_gaussian_fixed_mu(rt, mu_values, *popt_multi_gaussian),
            "r-",
            lw=1.5,
            label="Gaussian Fit",
        )
        ax2.set_xlim(0, 1.3)

    # Condition to plot if no peaks were detected or fitting fails
    if len(peak_indices) == 0 or popt_multi_gaussian == []:
        ax2.set_xlim(0, 1.3)
    else:
        """ax2.set_xlim(rt[peak_indices[0]] - 1, rt[peak_indices[-1]] + 1)"""
        ax2.set_xlim(0, 1.3)
    y_values = multi_gaussian_fixed_mu(
        rt[peak_indices], mu_values, *popt_multi_gaussian
    )

    # Add LC chromatogram peak annotations
    for j in peak_indices:
        label = str(rt[j])
        label_y_pos = max(rt_i) + 0.02 * max(rt_i)
        ax2.annotate(
            label[:4],
            xy=(rt[j], label_y_pos),
            fontsize=10,
            fontweight="bold",
            fontname="Arial",
            ha="center",
            va="bottom",
        )
    ax2.scatter(
        rt[peak_indices],
        y_values,
        color="purple",
        marker="*",
        s=40,
        label="Identified Peaks",
    )

    # Format LC chromatogram
    ax2.set_xlabel(
        "Retention Time [min]",
        fontsize=10,
        fontweight="bold",
        fontname="Arial",
    )
    ax2.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
    ax2.set_title(
        "Extracted Chromatogram",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
    )
    ax2.tick_params(axis="both", which="major", labelsize=10, width=1)
    for label in ax2.get_xticklabels():
        label.set_fontname("Arial")
        label.set_fontweight("bold")
        label.set_fontsize(10)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_linewidth(1.5)
    ax2.spines["bottom"].set_linewidth(1.5)
    legend2 = ax2.legend(
        loc="best",
        frameon=True,
        fontsize=10,
        edgecolor="black",
        facecolor="white",
    )
    for text in legend2.get_texts():
        text.set_fontname("Arial")
    max_intensity_y_limit = max(rt_i) + 0.1 * max(rt_i)
    y_tick_values = np.linspace(0, max_intensity_y_limit, num=10, endpoint=True)
    ax2.set_yticks(y_tick_values)
    tick_label_fontprops = {"weight": "bold", "family": "Arial", "size": 10}
    ax2.set_yticklabels([int(y) for y in y_tick_values], fontdict=tick_label_fontprops)
    ax2.set_ylim(0, max_intensity_y_limit)

    # Plot mobilogram
    peak_height = max(fit_i)
    peak_index = np.argmax(fit_i)
    peak_x_position = t_refined[peak_index]
    offset = peak_height * 0.001
    ax3.text(
        peak_x_position,
        peak_height + offset,
        "{:.2f}".format(B),
        c="k",
        fontsize=10,
        fontweight="bold",
        fontname="Arial",
        ha="center",
        va="bottom",
    )
    ax3.plot(t, dt_i, "o--b", lw=1.5, ms=2, label="Raw Data")
    ax3.plot(t_refined, fit_i, "r-", lw=1.5, label="Gaussian Fit")
    ax3.set_title(
        f"Extracted Mobilogram\nObserved m/z: {monoisotopic_mz:.4f} \u00B1 0.025 rt: {rt_value} → FWHM ~ {fwhm:.2f} ",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
    )
    ax3.set_xlabel(
        "Drift Time [ms]",
        fontsize=10,
        fontweight="bold",
        fontname="Arial",
    )
    ax3.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
    max_dt_i_y_limit = max(dt_i) + 0.1 * max(dt_i)
    ax3.set_ylim(0, max_dt_i_y_limit)
    ax3.set_xlim(0, B + 2)
    ax3.tick_params(axis="both", which="major", labelsize=10, width=1)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.spines["left"].set_linewidth(1.5)
    ax2.spines["bottom"].set_linewidth(1.5)
    legend3 = ax3.legend(
        loc="best",
        frameon=True,
        fontsize=10,
        edgecolor="black",
        facecolor="white",
    )
    for text in legend3.get_texts():
        text.set_fontname("Arial")

    # Generate chemical structure from SMILES string
    chem_struct_img = smiles_to_structure(smiles, img_size=(350, 350))

    # Convert PIL image to array so it can be displayed
    if chem_struct_img is not None:
        chem_struct_arr = np.array(chem_struct_img)

        # Insert the chemical structure above the title for ax1
        insert_ax = fig.add_axes([0.2, 0.7, 0.25, 0.25])
        insert_ax.imshow(chem_struct_arr)
        insert_ax.axis("off")

    # Plot MS1 scan
    ax1.plot(
        mz_spectrum,
        intensity_spectrum,
        lw=1.5,
        color="black",
        ms=2,
        label="Raw MS1 Data",
    )

    # Add annotations for monoisotopic peak
    max_y = monoisotopic_intensity * (1.1)
    ax1.text(
        monoisotopic_mz - 0.01,
        max_y * 0.98,
        f"{monoisotopic_mz:.4f}",
        fontsize=10,
        fontweight="bold",
        fontname="Arial",
        ha="right",
        va="top",
    )
    ax1.axvline(x=monoisotopic_mz, color="magenta", label="Monoisotopic Peak", lw=2.5)
    mz_min = monoisotopic_mz - 1
    mz_max = monoisotopic_mz + 1
    ax1.set_xlim(mz_min, mz_max)
    ax1.set_xlabel("m/z", fontsize=10, fontweight="bold", fontname="Arial")
    ax1.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
    ax1.set_title(
        f"{compound_name}\n{mz:.4f} {adduct}\nMS1 Spectrum",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
    )
    ax1.tick_params(axis="both", which="major", labelsize=10, width=1)
    y_ticks = np.linspace(0, max_y, num=10, endpoint=True)
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels([int(y) for y in y_ticks], fontdict=tick_label_fontprops)
    ax1.set_ylim(0, max_y)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_linewidth(1.5)
    ax1.spines["bottom"].set_linewidth(1.5)
    legend1 = ax1.legend(
        loc="upper right",
        frameon=True,
        fontsize=10,
        edgecolor="black",
        facecolor="white",
    )
    warnings.filterwarnings(
        "ignore",
        message="This figure includes Axes that are not compatible with tight_layout, so results might be incorrect",
    )
    for text in legend1.get_texts():
        text.set_fontname("Arial")

    # Attempt to save the figure
    try:
        plt.tight_layout()
        plt.savefig(fname_combined, dpi=300, bbox_inches="tight")
    except Exception as e:
        plt.savefig(fname_combined, dpi=300)
    plt.close()


def combined_multi_figure(
    rt,
    rt_i,
    t,
    dt_i,
    file_name,
    compound_name,
    adduct,
    mz,
    monoisotopic_mz,
    monoisotopic_intensity,
    rt_value,
    mu_values,
    fname_combined,
    mz_spectrum,
    intensity_spectrum,
    smiles,
):
    """
    toxccs/utils/plotting.py/combined_multi_figure
    description:
            Generates a figure containing the extracted chromatogram, mobilogram, and MS1 scan with the identified monoisotopic peak displayed. Specific logic included for handling complex ion mobilogram profiles.
    parameters:
            rt (list or numpy.ndarray) -- array of retention times.
            rt_i (list or numpy.ndarray) -- corresponding intensity values for each retention time.
            t (list or numpy.ndarray) -- raw EIM time points.
            dt_i (list or numpy.ndarray) -- raw EIM intensity values.
            t_refined (list or numpy.ndarray) -- refined data points for plotting the fit.
            fit_i (list or numpy.ndarray) -- intensity values of the fitted curve.
            A (float) -- Gaussian parameter.
            B (float) -- Gaussian parameter.
            file_name (str) --name of the .raw file from which the data is extracted.
            compound_name (str) -- name of the compound being analyzed.
            adduct (str) -- adduct type of the compound being analyzed.
            mz (float) -- target m/z value used to extract the raw data.
            monoisotopic_mz (float) -- observed monoisotopic m/z value.
            rt_value (float) -- observed retention time.
            fwhm (float) -- FWHM of the EIM.
            mu_values -- (ndarray) -- mu values for multi-Gaussian function.
            fname_combined (str) -- file name for saving the plot.
            smoothed_ms1_data (list or numpy.ndarray) -- smoothed MS1 data.
            smiles (str) -- SMILES string of the chemical being processed.
    returns:
            (.png) -- combined image of the extracted and fitted data.
    """

    window_len = 51
    window = gaussian(window_len, std=1)
    smoothed_rt_i = convolve(rt_i, window / window.sum(), mode="same")
    smoothed_dt_i = convolve(dt_i, window / window.sum(), mode="same")

    # Retrieve LC chromatogram information
    (
        valid_peak_indices,
        _,
        _,
        peak_ranges,
        popt_multi_gaussian,
        mu_values,
        _,
        _,
    ) = process_multi_chromatogram(
        rt,
        rt_i,
    )

    # Set up plot
    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(6.4, 9.6))

    # Plot LC chromatogram figure for successfully fitted EICs
    ax2.plot(rt, rt_i, "b-", lw=1.5, label="Raw Data")
    ax2.plot(rt, smoothed_rt_i, "g--", lw=1.5, label="Smoothed Data")

    # Condition to plot if peaks were successfully fitted
    if len(valid_peak_indices) > 0:
        y_values = multi_gaussian_fixed_mu(
            rt[valid_peak_indices], mu_values, *popt_multi_gaussian
        )
        ax2.plot(
            rt,
            multi_gaussian_fixed_mu(rt, mu_values, *popt_multi_gaussian),
            "r-",
            lw=1.5,
            label="Gaussian Fit",
        )

        # Visualize the retention time bounds as vertical lines
        """for start_idx, end_idx in peak_ranges:
                rt_start = rt[start_idx]
                rt_end = rt[end_idx]
                ax2.axvline(x=rt_start, color="blue", linestyle="--", linewidth=1.5)
                ax2.axvline(x=rt_end, color="red", linestyle="--", linewidth=1.5)

            # Plot dashed lines for FWHM for each peak
            for i in range(0, len(popt_multi_gaussian), 2):
                A = popt_multi_gaussian[i]  # Amplitude of the peak
                sigma = popt_multi_gaussian[i + 1]  # Standard deviation of the peak
                mu = mu_values[int(i / 2)]  # Mean (retention time apex)
                fwhm = sigma * 2.355
                half_max = A / 2.0
                lower_bound = mu - (fwhm / 2)
                upper_bound = mu + (fwhm / 2)

                # Plot horizontal line at half-max height, spanning the FWHM range
                ax2.hlines(
                    y=half_max,
                    xmin=lower_bound,
                    xmax=upper_bound,
                    colors="purple",
                    linestyles="solid",
                    linewidth=1.5,
                    label="FWHM",
                )

                # Annotate the FWHM value
                ax2.text(
                    mu,
                    half_max + 0.05 * max(rt_i),
                    f"FWHM: {fwhm:.3} min",
                    fontsize=8,
                    color="purple",
                    ha="center",
                    va="bottom",
                )"""

        # Fill the area under each Gaussian curve
        """for i in range(0, len(popt_multi_gaussian), 2):
                A = popt_multi_gaussian[i]
                sigma = popt_multi_gaussian[i + 1]
                mu = mu_values[int(i / 2)]

                print(
                    f"Peak {int(i/2) + 1} | Amplitude (A): {A}, Sigma: {sigma}, Mu: {mu}"
                )

                # Calculate the Gaussian curve for this peak
                gaussian_curve = A * np.exp(-((rt - mu) ** 2) / (2 * sigma**2))

                # Shade the area under the Gaussian curve
                ax2.fill_between(
                    rt,
                    gaussian_curve,
                    color="purple",
                    alpha=0.3,
                    label=f"Peak {int(i/2) + 1} Area",
                )"""

        ax2.set_xlim(0, 1.3)

    # Condition to plot if no peaks were detected or fitting fails
    if len(valid_peak_indices) == 0 or popt_multi_gaussian == []:
        ax2.set_xlim(0, 1.3)
        y_values = []
        valid_peak_indices = []
    else:
        ax2.set_xlim(0, 1.3)

    # Add LC chromatogram peak annotations
    for idx, peak_idx in enumerate(valid_peak_indices):
        label = str(rt[peak_idx])
        label_y_pos = y_values[idx] + 0.03 * max(y_values)
        """label_y_pos = rt_i[j] + 0.02 * max(rt_i)"""
        ax2.annotate(
            label[:4],
            xy=(rt[peak_idx], label_y_pos),
            fontsize=10,
            fontweight="bold",
            fontname="Arial",
            ha="center",
            va="bottom",
        )
    ax2.scatter(
        rt[valid_peak_indices],
        y_values,
        color="purple",
        marker="*",
        s=40,
        label="Identified Peaks",
    )

    # Format LC chromatogram
    ax2.set_xlabel(
        "Retention Time [min]",
        fontsize=10,
        fontweight="bold",
        fontname="Arial",
    )
    ax2.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
    ax2.set_title(
        "Extracted Chromatogram",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
    )
    ax2.tick_params(axis="both", which="major", labelsize=10, width=1)
    for label in ax2.get_xticklabels():
        label.set_fontname("Arial")
        label.set_fontweight("bold")
        label.set_fontsize(10)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.spines["left"].set_linewidth(1.5)
    ax2.spines["bottom"].set_linewidth(1.5)
    legend2 = ax2.legend(
        loc="best",
        frameon=True,
        fontsize=10,
        edgecolor="black",
        facecolor="white",
    )
    for text in legend2.get_texts():
        text.set_fontname("Arial")
    max_intensity_y_limit = max(rt_i) + 0.2 * max(rt_i)
    y_tick_values = np.linspace(0, max_intensity_y_limit, num=10, endpoint=True)
    ax2.set_yticks(y_tick_values)
    tick_label_fontprops = {"weight": "bold", "family": "Arial", "size": 10}
    ax2.set_yticklabels([int(y) for y in y_tick_values], fontdict=tick_label_fontprops)
    ax2.set_ylim(0, max_intensity_y_limit)

    # Retrieve mobilogram information
    (
        valid_dt_peak_indices,
        _,
        _,
        _,
        dt_popt_multi_gaussian,
        dt_mu_values,
        _,
        _,
    ) = process_multi_chromatogram(
        t,
        dt_i,
    )

    # Plot mobilogram figure for successfully fitted EIMs
    ax3.plot(t, dt_i, "b-", lw=1.5, label="Raw Data")
    ax3.plot(t, smoothed_dt_i, "g--", lw=1.5, label="Smoothed Data")
    ax3.set_title(
        f"Extracted Mobilogram\nObserved m/z: {monoisotopic_mz:.4f} \u00B1 0.025 rt: {rt_value}",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
    )

    # Define tmax for plotting
    max_dt_i = np.argmax(dt_i)
    max_dt = t[max_dt_i]

    # Condition to plot if peaks were successfully fitted
    if len(valid_dt_peak_indices) > 0:
        ax3.plot(
            t,
            multi_gaussian_fixed_mu(t, dt_mu_values, *dt_popt_multi_gaussian),
            "r-",
            lw=1.5,
            label="Gaussian Fit",
        )
        ax3.set_xlim(max_dt - 2, max_dt + 2)

    # Condition to plot if no peaks were detected or fitting fails
    if len(valid_dt_peak_indices) == 0 or dt_popt_multi_gaussian == []:
        ax3.set_xlim(max_dt - 2, max_dt + 2)
    else:
        ax3.set_xlim(max_dt - 2, max_dt + 2)
    dt_y_values = multi_gaussian_fixed_mu(
        t[valid_dt_peak_indices], dt_mu_values, *dt_popt_multi_gaussian
    )

    # Add mobilogram peak annotations
    for idx, peak_idx in enumerate(valid_dt_peak_indices):
        dt_label = str(t[peak_idx])
        label_dt_pos = dt_y_values[idx] + 0.03 * max(dt_y_values)
        ax3.annotate(
            dt_label[:4],
            xy=(t[peak_idx], label_dt_pos),
            fontsize=10,
            fontweight="bold",
            fontname="Arial",
            ha="center",
            va="bottom",
        )
    ax3.scatter(
        t[valid_dt_peak_indices],
        dt_y_values,
        color="purple",
        marker="*",
        s=40,
        label="Identified Peaks",
    )

    # Format mobilogram
    ax3.set_xlabel("Drift Time [ms]", fontsize=10, fontweight="bold", fontname="Arial")
    ax3.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
    ax3.tick_params(axis="both", which="major", labelsize=10, width=1)
    for label in ax3.get_xticklabels():
        label.set_fontname("Arial")
        label.set_fontweight("bold")
        label.set_fontsize(10)
    ax3.spines["top"].set_visible(False)
    ax3.spines["right"].set_visible(False)
    ax3.spines["left"].set_linewidth(1.5)
    ax3.spines["bottom"].set_linewidth(1.5)
    legend3 = ax3.legend(
        loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white"
    )
    for text in legend3.get_texts():
        text.set_fontname("Arial")
    max_intensity_y_limit_dt = max(dt_i) + 0.2 * max(dt_i)
    y_tick_values_dt = np.linspace(0, max_intensity_y_limit_dt, num=10, endpoint=True)
    ax3.set_yticks(y_tick_values_dt)
    tick_label_fontprops = {"weight": "bold", "family": "Arial", "size": 10}
    ax3.set_yticklabels(
        [int(y) for y in y_tick_values_dt], fontdict=tick_label_fontprops
    )
    ax3.set_ylim(0, max_intensity_y_limit_dt)

    # Generate chemical structure from SMILES string
    chem_struct_img = smiles_to_structure(smiles, img_size=(350, 350))

    # Convert PIL image to array so it can be displayed
    if chem_struct_img is not None:
        chem_struct_arr = np.array(chem_struct_img)

        # Insert the chemical structure above the title for ax1
        insert_ax = fig.add_axes([0.2, 0.7, 0.25, 0.25])
        insert_ax.imshow(chem_struct_arr)
        insert_ax.axis("off")

    # Plot MS1 scan
    ax1.plot(
        mz_spectrum,
        intensity_spectrum,
        lw=1.5,
        color="black",
        ms=2,
        label="Raw MS1 Data",
    )

    # Add annotations for monoisotopic peak
    max_y = monoisotopic_intensity * (1.1)
    ax1.text(
        monoisotopic_mz - 0.01,
        max_y * 0.98,
        f"{monoisotopic_mz:.4f}",
        fontsize=10,
        fontweight="bold",
        fontname="Arial",
        ha="right",
        va="top",
    )
    ax1.axvline(x=monoisotopic_mz, color="magenta", label="Monoisotopic Peak", lw=2.5)
    mz_min = monoisotopic_mz - 1
    mz_max = monoisotopic_mz + 1
    ax1.set_xlim(mz_min, mz_max)
    ax1.set_xlabel("m/z", fontsize=10, fontweight="bold", fontname="Arial")
    ax1.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
    ax1.set_title(
        f"{compound_name}\n{mz:.4f} {adduct}\nMS1 Spectrum",
        fontsize=12,
        fontweight="bold",
        fontname="Arial",
    )
    ax1.tick_params(axis="both", which="major", labelsize=10, width=1)
    y_ticks = np.linspace(0, max_y, num=10, endpoint=True)
    ax1.set_yticks(y_ticks)
    ax1.set_yticklabels([int(y) for y in y_ticks], fontdict=tick_label_fontprops)
    ax1.set_ylim(0, max_y)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_linewidth(1.5)
    ax1.spines["bottom"].set_linewidth(1.5)
    legend1 = ax1.legend(
        loc="upper right",
        frameon=True,
        fontsize=10,
        edgecolor="black",
        facecolor="white",
    )
    warnings.filterwarnings(
        "ignore",
        message="This figure includes Axes that are not compatible with tight_layout, so results might be incorrect",
    )
    for text in legend1.get_texts():
        text.set_fontname("Arial")

    # Attempt to save the figure
    try:
        plt.tight_layout()
        plt.savefig(fname_combined, dpi=300, bbox_inches="tight")
    except Exception as e:
        plt.savefig(fname_combined, dpi=300)
    plt.close()


def combined_dda_figure(
    rt,
    rt_i,
    ms2,
    ms2_i,
    file_name,
    compound_name,
    adduct,
    mz,
    fname_combined,
    smiles,
    rt_min,
    rt_max,
    x_min,
    x_max,
):
    """
    toxccs/utils/plotting.py/combined_dda_figure
    description:
            Generates a figure containing the extracted chromatogram and MS/MS fragmentation spectrum with the top 5 most intense m/z peaks highlighted.
    parameters:
            rt (list or numpy.ndarray) -- array of retention times.
            rt_i (list or numpy.ndarray) -- corresponding intensity values for each retention time.
            ms2 (list or numpy.ndarray) -- array of m/z values.
            ms2_i (list or numpy.ndarray) -- corresponding intensity values for each m/z.
            file_name (str) -- raw file name for creating the output file name.
            compound_name (str) -- name of the compound being analyzed.
            adduct (str) -- adduct type of the compound being analyzed.
            mz (float) -- target m/z value used to extract the raw data.
            fname_combined (str) -- file name for saving the plot.
            smiles (str) -- SMILES string of the chemical being processed.
            rt_min (float) -- lowest time bound of the extracted LC ion chromatogram.
            rt_max (float) -- highest time bound of the extracted LC ion chromatogram.
    returns:
            (.png) -- combined image of the extracted and fitted data.
    """
    try:

        # Ensure ms2 and ms2_i arrays are not empty
        if ms2.size == 0 or ms2_i.size == 0:
            raise ValueError("MS2 data arrays are empty.")

        # Smooth extracted LC ion chromatogram
        window_len = 51
        window = gaussian(window_len, std=1)
        smoothed_rt_i = convolve(rt_i, window / window.sum(), mode="same")

        # Retrieve LC chromatogram information
        peak_indices, rt_list, areas, peak_ranges, popt_multi_gaussian, mu_values = (
            process_chromatogram(
                rt,
                rt_i,
                file_name,
                compound_name,
                adduct,
                mz,
                "sample",
                generate_images=False,
            )
        )

        # Set up plot
        fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(6.4, 6.4))

        # Plot LC chromatogram figure for successfully fitted EICs
        ax1.plot(rt, rt_i, "b-", lw=1.5, label="Raw Data")
        ax1.plot(rt, smoothed_rt_i, "g--", lw=1.5, label="Smoothed Data")

        # Condition to plot if peaks were successfully fitted
        if len(peak_indices) > 0:
            ax1.plot(
                rt,
                multi_gaussian_fixed_mu(rt, mu_values, *popt_multi_gaussian),
                "r-",
                lw=1.5,
                label="Gaussian Fit",
            )
            ax1.set_xlim(x_min, x_max)

        # Condition to plot if no peaks were detected or fitting fails
        if len(peak_indices) == 0 or popt_multi_gaussian == []:
            ax1.set_xlim(x_min, x_max)
        else:
            ax1.set_xlim(x_min, x_max)
        y_values = multi_gaussian_fixed_mu(
            rt[peak_indices], mu_values, *popt_multi_gaussian
        )

        # Add LC chromatogram peak annotations
        for j in peak_indices:
            peak_height = rt_i[j]
            label = str(rt[j])
            label_y_pos = peak_height + 0.02 * peak_height
            # label_y_pos = max(rt_i) + 0.02 * max(rt_i)
            ax1.annotate(
                label[:4],
                xy=(rt[j], label_y_pos),
                fontsize=10,
                fontweight="bold",
                fontname="Arial",
                ha="center",
                va="bottom",
            )

        ax1.scatter(
            rt[peak_indices],
            rt_i[peak_indices],
            # y_values,
            color="purple",
            marker="*",
            s=40,
            label="Identified Peaks",
        )

        # Plot vertical lines for rt_min and rt_max time bounds
        ax1.axvline(x=rt_min, color="black", linestyle="--", linewidth=1.5)
        ax1.axvline(x=rt_max, color="black", linestyle="--", linewidth=1.5)

        # Format LC chromatogram
        ax1.set_xlabel(
            "Retention Time [min]", fontsize=10, fontweight="bold", fontname="Arial"
        )
        ax1.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
        ax1.set_title(
            f"{compound_name}\n{mz:.4f} {adduct}\nExtracted Chromatogram",
            fontsize=12,
            fontweight="bold",
            fontname="Arial",
        )
        ax1.tick_params(axis="both", which="major", labelsize=10, width=1)
        for label in ax1.get_xticklabels():
            label.set_fontname("Arial")
            label.set_fontweight("bold")
            label.set_fontsize(10)
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.spines["left"].set_linewidth(1.5)
        ax1.spines["bottom"].set_linewidth(1.5)
        legend1 = ax1.legend(
            loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white"
        )
        for text in legend1.get_texts():
            text.set_fontname("Arial")
        max_intensity_y_limit = max(rt_i) + 0.1 * max(rt_i)
        y_tick_values = np.linspace(0, max_intensity_y_limit, num=10, endpoint=True)
        ax1.set_yticks(y_tick_values)
        tick_label_fontprops = {"weight": "bold", "family": "Arial", "size": 10}
        fontdict = tick_label_fontprops
        ax1.set_ylim(0, max_intensity_y_limit)

        # Generate chemical structure from SMILES string
        chem_struct_img = smiles_to_structure(smiles, img_size=(350, 350))

        # Convert PIL image to array so it can be displayed
        if chem_struct_img is not None:
            chem_struct_arr = np.array(chem_struct_img)

            # Insert the chemical structure above the title for ax1
            insert_ax = fig.add_axes([0.2, 0.65, 0.23, 0.23])
            insert_ax.imshow(chem_struct_arr)
            insert_ax.axis("off")

        # Plot extracted MS/MS fragmentation spectrum
        ax2.plot(ms2, ms2_i, lw=1.5, color="black", ms=2, label="Raw Data")
        mz_min = 0
        mz_max = mz + 10
        ax2.set_xlim(mz_min, mz_max)

        # Highlight top 5 most intense peaks
        in_range_mask = (ms2 >= mz_min) & (ms2 <= mz_max)
        ms2_in_range = ms2[in_range_mask]
        ms2_i_in_range = ms2_i[in_range_mask]
        top_5_indices = np.argsort(ms2_i_in_range)[-5:][::-1]
        top_5_mz = ms2_in_range[top_5_indices]
        top_5_intensities = ms2_i_in_range[top_5_indices]
        for peak_mz, peak_intensity in zip(top_5_mz, top_5_intensities):
            ax2.axvline(
                x=peak_mz,
                ymin=0,
                ymax=peak_intensity / ms2_i_in_range.max(),
                color="magenta",
                lw=2,
            )
            max_y = peak_intensity * (1.1)
            x_offset = 0.005 * (mz_max - mz_min)
            ax2.text(
                peak_mz - x_offset,
                max_y * 0.91,
                f"{peak_mz:.4f}",
                fontsize=8,
                fontweight="bold",
                fontname="Arial",
                ha="right",
                va="top",
            )
        ax2.axvline(x=0, color="magenta", lw=2, label="Top 5 Peaks")
        ax2.set_xlabel("m/z", fontsize=10, fontweight="bold", fontname="Arial")
        ax2.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
        ax2.set_title(
            f"Extracted MS/MS Fragmentation Specrum",
            fontsize=12,
            fontweight="bold",
            fontname="Arial",
        )
        ax2.tick_params(axis="both", which="major", labelsize=10, width=1)
        y_ticks = np.linspace(0, ms2_i.max(), num=10, endpoint=True)
        ax2.set_yticks(y_ticks)
        ax2.set_yticklabels([int(y) for y in y_ticks], fontdict=tick_label_fontprops)
        ax2.set_ylim(0, ms2_i.max())
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.spines["left"].set_linewidth(1.5)
        ax2.spines["bottom"].set_linewidth(1.5)
        legend2 = ax2.legend(
            loc="best",
            frameon=True,
            fontsize=10,
            edgecolor="black",
            facecolor="white",
        )
        warnings.filterwarnings(
            "ignore",
            message="This figure includes Axes that are not compatible with tight_layout, so results might be incorrect",
        )
        for text in legend2.get_texts():
            text.set_fontname("Arial")

        # Attempt to save the figure
        try:
            plt.tight_layout()
            plt.savefig(fname_combined, dpi=300, bbox_inches="tight")
        except Exception as e:
            fig.set_size_inches(10, 8)
            plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, hspace=0.5)
            plt.savefig(fname_combined, dpi=300)
        plt.close()

    except Exception as e:
        print(f"Error generating figure for {compound_name} ({file_name}): {e}")
