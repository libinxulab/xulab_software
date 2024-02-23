"""
    toxccs/extract.py
    Ryan Nguyen
    2/20/2024

    description:
    Module designed to handle the extraction and basic processing of full-scan LC-IM-MS/MS data. Provides a comprehensive set of tools for extracting raw data from .raw files, detecting peaks in the time and m/z dimensions, and converting drift times to calibrated collision cross-section (CCS) values. 
"""

# Todo
# Test monoisotopic peak algorithm
# Modify extract logic so that observed m/z is used to extract rt
# Add MS1 level scan to combined figure

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, gaussian, convolve
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from dhrmasslynxapi.reader import MassLynxReader
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRawXL

# Set global font conditions for figures
params = {"font.size": 8, "font.family": "Arial", "font.weight": "bold"}
plt.rcParams.update(params)


class RawProcessor:
    def __init__(self, target_list, reference_is_file=None):
        """
        RawProcessor.__init__
        description:
                Initializes a new RawProcessor with a target list containing precursor ions and
                paths to .raw data files.
        parameters:
                target_list (str) -- path to the Excel file (.xlsx) containing the target list.
                reference_is_list (str) -- path to th Excel file (.xlsx) containing the reference internal standard list (optional).
        """
        self.target_list = target_list
        self.reference_is_df = None
        if reference_is_file:
            self.read_reference_is(reference_is_file)
        self.output_rows = []

    def multi_gaussian_fixed_mu(self, x, mu_values, *params):
        """
        RawProcessor.multi_gaussian_fixed_mu
        description:
                Computes the sum of multiple Gaussian functions with fixed mean (mu) values for each peak. This function is primarily used for fitting chromatographic data where the mean positions of the peaks are known, but their ampltudes and widths (standard deviations) are to be determined.
        parameters:
                x (numpy.ndarray) -- independent variable (e.g., rt or m/z) where the function is evaluated.
                mu_values (numpy.ndarray) -- array of fixed mean (mu) values for the Gaussian functions. Each value corresponds to the mean position of a peak in the data.
                *params (list) -- parameters for the Gaussian functions, with each peak's amplitude and standard deviation.
        returns:
                y (numpy.ndarray) -- the evaluated sum of Gaussian functions at each point in x.
        """
        y = np.zeros_like(x)
        for i in range(0, len(params), 2):
            A = params[i]
            sigma = params[i + 1]
            mu = mu_values[int(i / 2)] if len(mu_values) > int(i / 2) else 0
            y += A * np.exp(-((x - mu) ** 2) / (2 * sigma**2))
        return y

    def process_chromatogram(
        self,
        rt,
        rt_i,
        file_name,
        compound_name,
        adduct,
        mz,
        sample_type,
        generate_images=False,
        prominence=1000,
        distance=3,
    ):
        """
        RawProcessor.process_chromatogram
        description:
            Processes raw data by smoothing time or m/z-intensity pairs using Gaussian convolution, identifying peaks, and
            fitting the smoothed data to a multi-Gaussian function. Optionally generates and saves chromatogram images.
        parameters:
            rt (list or numpy.ndarray) -- array of retention times.
            rt_i (list or numpy.ndarray) -- corresponding intensity values for each retention time.
            file_name (str) -- name of the .raw file from which the data is extracted.
            mz (float) -- target m/z value used to extract the raw data.
            sample_type (str) -- sample type (i.e., "control," "sample," or "IS")
            generation_images (bool, optional) -- flag to determine whether to generate and save chromatogram images.
            Default is false.
        returns:
            (tuple of (list of int, list of float, list of float, list of tuples))
                -- peak_indices -- indices of the identified peaks in the smoothed data.
                -- rt_list -- list of peaks corresponding to the identified peaks, rounded to 2 decimal places.
                -- areas -- calculated areas under each peak based on the multi-Gaussian fitting.
                -- peak_ranges -- start and end indices for each identified peak
        """

        window_len = 51
        window = gaussian(window_len, std=1)  # Adjust to change degree of smoothing
        smoothed_intensity = convolve(rt_i, window / window.sum(), mode="same")

        # Identify peaks in the smoothed data
        peak_indices, _ = find_peaks(
            smoothed_intensity, prominence=prominence, distance=distance
        )  # Adjust peak picking filters (default: prominence=1000, distance=3)

        # Identify peak indices
        peak_ranges = []
        for peak_idx in peak_indices:

            # Find the start and end of the peak
            start_idx = peak_idx
            while start_idx > 0 and rt_i[start_idx - 1] < rt_i[start_idx]:
                start_idx -= 1
            end_idx = peak_idx
            while end_idx < len(rt_i) - 1 and rt_i[end_idx + 1] < rt_i[end_idx]:
                end_idx += 1

            # Handle cases where start_idx is equal to or greater than peak apex
            if start_idx >= peak_idx:

                # Scan earlier points to find the appropriate start index
                new_start_idx = start_idx - 1
                while (
                    new_start_idx > 0 and rt_i[new_start_idx - 1] <= rt_i[new_start_idx]
                ):
                    new_start_idx -= 1
                start_idx = new_start_idx

            # Handle cases where end_idx is equal to or less than peak apex
            if end_idx <= peak_idx:

                # Scan later points to find the appropriate end index
                new_end_idx = end_idx + 1
                while (
                    new_end_idx < len(rt_i) - 1
                    and rt_i[new_end_idx + 1] <= rt_i[new_end_idx]
                ):
                    new_end_idx += 1
                end_idx = new_end_idx
            peak_ranges.append((start_idx, end_idx))

        # Add identified peaks to rt_list for output file
        rt_list = []
        for j in peak_indices:
            label = float(rt[j])
            round_label = str(round(label, 2))
            rt_list.append(round_label)

        # Extract mu (mean) values for the fixed Gaussian functions from the smoothed data
        mu_values = rt[peak_indices] if peak_indices.size > 0 else np.array([])

        # Initial guesses for amplitude (A) and sigma from the smoothed data
        A_guesses = (
            smoothed_intensity[peak_indices] if peak_indices.size > 0 else np.array([])
        )
        sigma_guesses = [0.05] * len(peak_indices)
        p0 = [val for sublist in zip(A_guesses, sigma_guesses) for val in sublist]

        # Attempt to fit the smoothed data to a multi-Gaussian function
        areas = []
        try:
            if peak_indices.size > 0:
                popt_multi_gaussian, _ = curve_fit(
                    lambda x, *p: self.multi_gaussian_fixed_mu(x, mu_values, *p),
                    rt,
                    rt_i,
                    p0=p0,
                    maxfev=20000,
                )
                for i in range(0, len(popt_multi_gaussian), 2):
                    A = popt_multi_gaussian[i]
                    sigma = popt_multi_gaussian[i + 1]
                    mu = mu_values[int(i / 2)]
                    area = A * sigma * np.sqrt(2 * np.pi)  # Gaussian integral
                    areas.append(area)
            else:
                popt_multi_gaussian = []
        except Exception as e:
            areas = [None] * len(peak_indices)
            popt_multi_gaussian = []

        return peak_indices, rt_list, areas, peak_ranges, popt_multi_gaussian, mu_values

    def peak_fit(self, t, dt_i, p0="guess"):
        """
        RawProcessor.peak_fit
        description:
                Fits a peak in the EIM to a Gaussian curve.
        parameters:
                t (list or numpy.ndarray) -- drift times.
                dt_i (list or numpy.ndarray) -- intensity  values corresponding to the drift times.
                p0 (tuple/str) -- initial guesses for the Gaussian parameters.
        returns:
                (array) -- optimal values for the Gaussian parameters.
        """
        if p0 == "guess":
            p0 = (max(dt_i), t[np.argmax(dt_i)], 0.5)
        opt, cov = curve_fit(self.gaussian_fit, t, dt_i, maxfev=5000, p0=p0)

        return opt

    def gaussian_fit(self, x, A, B, C):
        """
        RawProcessor.gaussian_fit
        description:
                Gaussian function used for fitting EIMs.
        parameters:
                x (array) -- data points where the Gaussian function is evaluated.
                A (float) -- amplitude of the Gaussian curve.
                B (float) -- mean value of the Gaussian curve.
                C (float) -- standard deviation of the Gaussian curve.
        returns:
                (array) -- evaluated Gaussian function at each point in x.
        """
        if abs(C) < 0.01:
            C = 0.01

        return A * np.exp(-((x - B) ** 2) / (2 * C**2))

    def fwhm_threshold(
        self, C, dt_i, fwhm_thresholds=(0.05, 2.5), intensity_threshold=500
    ):
        """
        RawProcessor.fwhm_threshold
        description:
                Checks if the full width at half maximum (FWHM) and intensity of the EIM meet
                threshold values.
        parameters:
                C (float) -- standard deviation of the Gaussian curve.
                dt_i (array) -- intensity values.
                fwhm_thresholds (tuple) -- minimum and maximum acceptable FWHM.
                intensity_threshold (float) -- minimum acceptale intensity.
        returns:
                (bool) -- True if the FWHM and intensity meet the thresholds, False otherwise.
        """
        fwhm = C * 2.355
        max_intensity = max(dt_i)

        return (
            fwhm_thresholds[0] < fwhm < fwhm_thresholds[1]
            and max_intensity > intensity_threshold
        )

    def atd(self, t, t_refined, dt_i, fit_i, A, B, title_atd, fname_atd):
        """
        RawProcessor.atd
        description:
                Generates a plot of the raw EIM and the Gaussian fit.
        parameters:
                t (list or numpy.ndarray) -- raw EIM time points.
                t_refined (list or numpy.ndarray) -- refined data points for plotting the fit.
                dt_i (list or numpy.ndarray) -- original intensity values.
                fit_i (list or numpy.ndarray) -- intensity values of the fitted curve.
                A, B (floats) -- Gaussian parameters.
                title_atd (str) -- title for the plot.
                fname_atd (str) -- file name for saving the plot.
        """
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.plot(t, dt_i, "o--b", lw=1.5, ms=2, label="Raw Data")  # Plot raw EIM
        ax.plot(
            t_refined, fit_i, "r-", lw=1.5, label="Gaussian Fit"
        )  # Plot fitted data
        legend = ax.legend(
            loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white"
        )
        for text in legend.get_texts():
            text.set_fontname("Arial")
        ax.text(
            B + 0.1,
            0.95 * max(fit_i),
            "{:.2f}".format(B),
            c="k",
            fontsize=10,
            fontweight="bold",
            fontname="Arial",
        )  # Annotate peak of fitted data with dt
        ax.set_title(title_atd, fontsize=12, fontweight="bold", fontname="Arial")
        ax.set_xlabel(
            "Drift Time [ms]", fontsize=12, fontweight="bold", fontname="Arial"
        )
        ax.set_ylabel("Intensity", fontsize=12, fontweight="bold", fontname="Arial")
        max_dt_i_y_limit = max(dt_i) + 0.1 * max(dt_i)
        plt.ylim(0, max_dt_i_y_limit)  # Set y axis range
        plt.xlim(0, B + 2)  # Set x axis range
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
        plt.close

    def combined_figure(
        self,
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
        rt_value,
        fwhm,
        mu_values,
        fname_combined,
    ):
        """
        RawProcessor.combined_figure
        description:
                Generates a figure containing the extracted chromatogram and extracted mobilogram.
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
        """
        window_len = 51
        window = gaussian(window_len, std=1)
        smoothed_rt_i = convolve(rt_i, window / window.sum(), mode="same")

        # Retrieve LC chromatogram information
        peak_indices, rt_list, areas, peak_ranges, popt_multi_gaussian, mu_values = (
            self.process_chromatogram(
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
                self.multi_gaussian_fixed_mu(rt, mu_values, *popt_multi_gaussian),
                "r-",
                lw=1.5,
                label="Gaussian Fit",
            )

        # Condition to plot if no peaks were detected or fitting fails
        if len(peak_indices) == 0 or popt_multi_gaussian == []:
            max_rt = rt[np.argmax(rt_i)]
            ax1.set_xlim(0, max_rt + 1)
        else:
            ax1.set_xlim(rt[peak_indices[0]] - 1, rt[peak_indices[-1]] + 1)
        y_values = self.multi_gaussian_fixed_mu(
            rt[peak_indices], mu_values, *popt_multi_gaussian
        )

        # Add LC chromatogram peak annotations
        for j in peak_indices:
            label = str(rt[j])
            ax1.annotate(
                label[:4],
                xy=(rt[j], rt_i[j]),
                xytext=(rt[j] + 0.04, rt_i[j] * 0.95),
                fontsize=10,
                fontweight="bold",
                fontname="Arial",
            )
        ax1.scatter(
            rt[peak_indices],
            y_values,
            color="purple",
            marker="*",
            s=40,
            label="Identified Peaks",
        )

        # Format LC chromatogram
        ax1.set_xlabel(
            "Retention Time [min]",
            fontsize=10,
            fontweight="bold",
            fontname="Arial",
        )
        ax1.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
        ax1.set_title(
            f"{compound_name} \nm/z: {mz:.4f} {adduct}\nExtracted Chromatogram",
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
        ax1.spines["right"].set_linewidth(1.5)
        ax1.spines["bottom"].set_linewidth(1.5)
        legend1 = ax1.legend(
            loc="best",
            frameon=True,
            fontsize=10,
            edgecolor="black",
            facecolor="white",
        )
        for text in legend1.get_texts():
            text.set_fontname("Arial")
        max_intensity_y_limit = max(rt_i) + 0.1 * max(rt_i)
        y_tick_values = np.linspace(0, max_intensity_y_limit, num=10, endpoint=True)
        ax1.set_yticks(y_tick_values)
        tick_label_fontprops = {"weight": "bold", "family": "Arial", "size": 10}
        ax1.set_yticklabels(
            [int(y) for y in y_tick_values], fontdict=tick_label_fontprops
        )
        ax1.set_ylim(0, max_intensity_y_limit)

        # Plot mobilogram
        ax2.text(
            B + 0.1,
            0.95 * max(fit_i),
            "{:.2f}".format(B),
            c="k",
            fontsize=10,
            fontweight="bold",
            fontname="Arial",
        )
        ax2.plot(t, dt_i, "o--b", lw=1.5, ms=2, label="Raw Data")
        ax2.plot(t_refined, fit_i, "r-", lw=1.5, label="Gaussian Fit")
        # Need to edit this still!
        ax2.set_title(
            f"Extracted Mobilogram\nObserved m/z: {monoisotopic_mz:.4f} \u00B1 0.025 rt: {rt_value} → FWHM ~ {fwhm:.2f} ",
            fontsize=12,
            fontweight="bold",
            fontname="Arial",
        )
        ax2.set_xlabel(
            "Drift Time [ms]",
            fontsize=10,
            fontweight="bold",
            fontname="Arial",
        )
        ax2.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
        max_dt_i_y_limit = max(dt_i) + 0.1 * max(dt_i)
        ax2.set_ylim(0, max_dt_i_y_limit)
        ax2.set_xlim(0, B + 2)
        ax2.tick_params(axis="both", which="major", labelsize=10, width=1)
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
        plt.tight_layout()
        plt.savefig(fname_combined, dpi=300, bbox_inches="tight")
        plt.close()

    def export_to_excel(self):
        """
        RawProcessor.export_to_excel
        description:
                Exports the extracted data to an Excel file (.xlsx).
        returns:
        (str) -- path to the output Excel file (.xlsx).
        """
        df_output = pd.DataFrame(
            self.output_rows,
            columns=[
                "File Name",
                "Compound Name",
                "Adduct",
                "Target m/z",
                "Observed m/z",
                "Column Type",
                "Sample Type" "CCS Calibrant",
                "Gradient",
                "Column Type",
                "Observed Drift Time (ms)",
                "Observed CCS (Å²)",
                "Observed Retention Time (min)",
                "EIC Peak Intensity",
                "EIC Peak Area",
            ],
        )
        base_name = os.path.splitext(os.path.basename(self.target_list))[0]
        output_file = "{}_processed.xlsx".format(base_name)
        df_output.to_excel(output_file, index=False)
        print("\nProcessing successful. View results in {}.".format(output_file))

        return output_file

    def read_reference_is(self, reference_is_file):
        """
        RawProcessor.read_reference_is
        description:
                Reads and stores the reference internal standards Excel file.
        parameters:
                reference_file (str) -- path to the Excel file (.xlsx) containing the reference IS list.
        """
        self.reference_is_df = pd.read_excel(reference_is_file)
        self.reference_is_dict = {}
        for _, row in self.reference_is_df.iterrows():
            key = (row["Exact m/z"], row["Gradient"])
            self.reference_is_dict[key] = {
                "Reference Retention Time (min)": row["Reference Retention Time (min)"],
                "Reference CCS (Å²)": row["Reference CCS (Å²)"],
            }

    def filter_data(self):
        """
        RawProcessor.filter_data
        description:
                Filters the extracted data to identify unique spectral features based on LC peak intensity.
        returns:
        (str) -- The path to the filtered output Excel file (.xlsx).
        """
        # Read the input Excel file containing extracted data, i.e., the output of RawProcess.extract
        df = pd.read_excel(self.output_file)

        # Convert retention time values to floats and strip out any flags
        df["Observed Retention Time (min)"] = df["Observed Retention Time (min)"].apply(
            lambda x: float(x.split()[0]) if isinstance(x, str) else x
        )

        # Select rows where the sample type is "sample"
        sample_rows = df[df["Sample Type"] == "sample"]

        # Initialize empty list to store rows that meet the filtering criteria
        unique_features = []

        # Select rows where the same type is "control" for comparison
        control_rows = df[df["Sample Type"] == "control"]

        # Iterate over each "sample" row to determine if it meets the filtering criteria
        for index, sample_row in sample_rows.iterrows():
            sample_mz = sample_row["Target m/z"]
            sample_rt = sample_row["Observed Retention Time (min)"]
            sample_peak_area = sample_row["EIC Peak Area"]

            # Start with an empty DataFrame for matching control features
            matching_controls = pd.DataFrame()

            # Check if there are any control features with the same m/z
            matching_mz_controls = control_rows[control_rows["Target m/z"] == sample_mz]

            # If there are matching m/z control features, further filter them by retention time
            if not matching_mz_controls.empty:
                matching_controls = matching_mz_controls[
                    abs(
                        matching_mz_controls["Observed Retention Time (min)"]
                        - sample_rt
                    )
                    <= 0.1
                ]

            # Handle case where there are no control features that match by m/z and retention time, or if the sample peak area is greater than the max control feature peak
            if matching_controls.empty or (
                sample_peak_area > 1.5 * matching_controls["EIC Peak Area"].max()
            ):

                # Assign a peak area ratio if matching control features are found
                if not matching_controls.empty:
                    sample_row["Control EIC Peak Area Ratio"] = (
                        sample_peak_area / matching_controls["EIC Peak Area"].max()
                    )
                else:

                    # Handle case where there are no matching control features
                    sample_row["Control EIC Peak Area Ratio"] = ""

                # Add the sample row to the list of unique features
                unique_features.append(sample_row)

        # Create a pandas DataFrame with the unique features
        df_output = pd.DataFrame(unique_features)

        # Define the output file name
        filtered_file = self.output_file.replace("_processed.xlsx", "_filtered.xlsx")

        # Export the DataFrame to an Excel file
        df_output.to_excel(filtered_file, index=False)
        print("\nFiltering successful. View results in {}\n".format(filtered_file))

        return filtered_file

    def gaussian_smooth_pick(
        self, mz_array, intensity_array, window_len=1, std=0.1, prominence=0.001
    ):
        """
        FeatureAnnotate.gaussian_smooth_pick
        description:
                Applies a Gaussian smoothing to the intensity array and identifies peaks.
                This function convolves the intensity array with a Gaussian window to smooth the data.
                It then identifies local maxima (i.e., peaks) in the smoothed data based on a defined prominence.
        parameters:
                mz_array (array-like): array of m/z values corresponding to the intensities.
                intensity_array (array-like): array of intensity values to be smoothed and from which peaks are identified.
                window_len (int) -- length of the gaussian window used for smoothing. Default is 1.
                std (float) -- standard deviation of the Gaussian window, controlling the degree of smoothing. Default is 0.1.
                prominence (float) -- minimum prominence of peaks to be identified. Default is 0.001
        returns:
                (tuple) -- a tuple containing identified peaks and smoothed intensity array.
        """
        window = gaussian(window_len, std=std)
        smoothed_intensity = convolve(
            intensity_array, window / window.sum(), mode="same"
        )
        peaks, _ = find_peaks(smoothed_intensity, prominence=prominence)
        identified_peaks = [(mz_array[p], smoothed_intensity[p]) for p in peaks]

        return identified_peaks, smoothed_intensity

    def monoisotopic_peak(self, identified_peaks, target_mz, tolerance=0.025):
        """
        RawProcessor.monoisotopic_peak
        description:
                Identifies the monoisotopic peak within the target rt-selected MS1 spectrum.
                The current implementation only considers the M+0 peak because all [M]+ analytes carry a permanent positive charge.
        parameters:
                identified_peaks (tuple) -- extracted and "centroided" MS1 peaks
                target_mz (float) -- target precursor m/z value to compare with extracted peaks.
                tolerance (float) -- tolerance around the target m/z value to consider.
        returns:
                (float) -- extracted m/z value of the monoisotopic peak.
        """
        # Identify the peaks from the extracted MS1 function that fall within the M+3 window
        potential_peaks = [
            (mz, intensity)
            for mz, intensity in identified_peaks
            if target_mz - tolerance <= mz <= target_mz + tolerance
        ]
        if not potential_peaks:
            return target_mz

        # Find the peak with the highest intensity within the tolerance window
        highest_intensity_peak = max(potential_peaks, key=lambda x: x[1])

        return highest_intensity_peak[0]

    def observed_mz(self, identified_peaks, target_mz, mz_tolerance=0.025):
        """
        RawProcessor.observed_mz
        description:
                A generalized version of RawProcessor.monoisotopic_mz. Identifies the observed (i.e., feature) m/z peak and isotopologue type within the MS1 scan. This function iteratively compares the intensity of the target feature with peaks at m/z - 1 Da, m/z - 2 Da, etc. The peak with the highest intensity is assumed to be the monoisotopic peak.
        parameters:
                identified_peaks  (list of tuples) -- extracted and "centroided" MS1 peaks
                target_mz (float) -- target precursor m/z value.
                mz_tolerance (float) -- tolerance around the target m/z value for identifying potential peaks.
        returns:
                (float) -- extracted m/z value of the monoisotopic peak.
        """

        # Filter peaks within the m/z tolerance range
        potential_peaks = [
            peak
            for peak in identified_peaks
            if abs(peak[0] - target_mz) <= mz_tolerance
        ]

        # If no peaks are found within tolerance, return the target m/z
        if not potential_peaks:
            return target_mz, 0

        # Assume the most intense peak is the monoisotopic peak
        monoisotopic_peak = max(potential_peaks, key=lambda x: x[1])
        print(monoisotopic_peak)

        # Initialize variables to hold the monoisotopic peak and its isotopologue type
        isotopologue_type = 0

        # Check for isotopologue peaks(M + 1, M + 2, ...) and compare intensities
        # Check up to M + 3 isotopologues
        for delta_mz in range(1, 4):
            isotopologue_mz = target_mz - delta_mz
            isotopologue_peak = next(
                (
                    peak
                    for peak in potential_peaks
                    if abs(peak[0] - isotopologue_mz) <= mz_tolerance
                ),
                None,
            )

            # Update the monoisotopic peak if an isotopologue is more intense
            if isotopologue_peak and isotopologue_peak[1] > monoisotopic_peak[1]:
                monoisotopic_peak = isotopologue_peak
                isotopologue_type = delta_mz

        return monoisotopic_peak[0], isotopologue_type

    def extract(self, calibration_file, ms1_function, mobility_function, mz_tolerance):
        """
        RawProcessor.etract
        description:
                Performs the extraction and processing of raw LC-IM-MS data.
        parameters:
                calibration_file (str) -- path to the CCS calibration file.
                ms1_function (int) -- MS1 function number.
                mobility_function (int) -- mobility function number.
                mz_tolerance (float) -- m/z tolerance value to extract data.
        """
        # Set the input parameters
        self.calibration_file = calibration_file
        self.ms1_function = ms1_function
        self.mobility_function = mobility_function
        self.mz_tolerance = mz_tolerance
        self.cal_data = CCSCalibrationRawXL(self.calibration_file)

        # Read target list into a pandas DataFrame
        df_input = pd.read_excel(self.target_list)
        df_input.reset_index(drop=True, inplace=True)

        # Generate a new folder called "Extracted Mobilograms" in the directory if not already present
        """mobilogram_directory = "Extracted Mobilograms"
        if not os.path.exists(mobilogram_directory):
            os.makedirs(mobilogram_directory)"""

        # Generate a new folder called "Extracted Data" in the directory if not already present
        data_directory = "Extracted Data"
        if not os.path.exists(data_directory):
            os.makedirs(data_directory)

        # Initialize dictionary for internal standard CCS flags
        ccs_flags = {}

        # Iterate over each row in the target list DataFrame
        print(
            "\n...CCS Calibration successful. Extracting data from the raw files using precursor target list {}...".format(
                self.target_list
            )
        )
        for i, row in df_input.iterrows():

            # Extract the path to .raw file, m/z, and sample type (i.e., sample or blank)
            file_name, mz, sample_type, gradient, adduct, compound_name = (
                row["File Name"],
                row["Target m/z"],
                row["Sample Type"],
                row["Gradient"],
                row["Adduct"],
                row["Compound Name"],
            )

            # Print the current m/z being queried to the terminal
            print(
                "\n\traw file: {} compound: {} m/z: {}".format(
                    file_name, compound_name, mz
                )
            )

            # Identify the monoisotopic peak within the MS1 scan
            # Average all scans within the acquisition time
            # Default is 0.2 to 1.15 min
            rdr = MassLynxReader(file_name)
            mz_spectrum, intensity_spectrum = rdr.get_spectrum(ms1_function, 0.2, 1.15)

            # Smooth and fit the extracted MS1 spectrum using Gaussian convolution
            # This process is analogous to "centroiding" the profile data as outlined in MSnbase
            identified_peaks, smoothed_intensity = self.gaussian_smooth_pick(
                mz_spectrum, intensity_spectrum
            )

            # Identify the monoisotopic peak in the MS1 centroided data
            monoisotopic_mz, _ = self.observed_mz(
                identified_peaks, mz, self.mz_tolerance
            )

            print(f"monoisotopic peak: {monoisotopic_mz}")

            # Extract m/z-selected ion chromatogram (EIC) from the .raw file using theoretical m/z value
            # Requires user-inputted MS1 function number and desired m/z tolerance
            rdr = MassLynxReader(file_name)
            rt, rt_i = rdr.get_chrom(
                self.ms1_function, float(monoisotopic_mz), self.mz_tolerance
            )

            # Store extracted rt, rt_i data as arrays for further processing
            rt = np.array(rt)
            rt_i = np.array(rt_i)

            # Smooth and fit EIC with multigauss module
            (
                peak_indices,
                rt_list,
                areas,
                peak_ranges,
                popt_multi_gaussian,
                mu_values,
            ) = self.process_chromatogram(
                rt,
                rt_i,
                file_name,
                compound_name,
                adduct,
                mz,
                sample_type,
                generate_images=True,
            )

            # Iterate over each extracted LC ion chromatogram peak
            for rt_value, peak_area, (start_idx, end_idx) in zip(
                rt_list, areas, peak_ranges
            ):
                if sample_type == "IS":

                    # Create a key for looking up the reference retention time based on m/z and gradient type
                    ref_rt_key = (mz, gradient)

                    # If the key is in the reference dictionary, check against the reference retention time
                    if ref_rt_key in self.reference_is_dict:
                        reference_rt = self.reference_is_dict[ref_rt_key][
                            "Reference Retention Time (min)"
                        ]

                        # Calculate the difference between observed and reference retention times
                        rt_difference = abs(float(rt_value) - reference_rt)

                        # Flag the retention time if difference is greater than 0.2 min
                        if rt_difference > 0.2:
                            rt_value = f"{rt_value} FLAG"

                # Use peak indices to define the window
                rt_start = rt[start_idx]
                rt_end = rt[end_idx]

                # Handle cases where peak indices are determined to be equal due to close spacing of raw data
                # Default to use a fixed 0.3 min window (i.e., +/- 0.15 min around the apex)
                if rt_start == rt_end:
                    fixed_bound = 0.15

                    # Adjust rt_start and rt_end for these cases
                    rt_start = max(rt[0], rt[start_idx] - fixed_bound)
                    rt_end = min(rt[-1], rt[end_idx] + fixed_bound)

                """# Extract MS1 spectrum for isotopologue check
                mz_spectrum, intensity_spectrum = rdr.get_spectrum(
                    ms1_function, rt_start, rt_end
                )

                # Smooth and fit the extracted MS1 spectrum using Gaussian convolution
                # This process is analogous to "centroiding" the profile data as outlined in MSnbase
                # The default parmaeters for the smoothing and picking (i.e., window_len=1, std=7, prominence=0.01) select the most intense peak for each ion distribution
                # More aggressive smoothing will lead to slightly different m/z peaks being picked
                identified_peaks, smoothed_intensity = self.gaussian_smooth_pick(
                    mz_spectrum, intensity_spectrum
                )

                # Identify the monoisotopic peak in the MS1 centroided data
                monoisotopic_mz = self.monoisotopic_peak(identified_peaks, float(mz))"""

                # Extract (m/z,rt)-selected ion mobilogram (EIM) for each identified LC peak
                t, dt_i = rdr.get_filtered_chrom(
                    self.mobility_function,
                    float(monoisotopic_mz),
                    self.mz_tolerance,
                    rt_min=rt_start,
                    rt_max=rt_end,
                )
                dt = None

                # Smooth and fit EIM to Gaussian function
                A, B, C = self.peak_fit(t, dt_i)
                t_refined = np.arange(min(t), max(t), 0.01)
                fit_i = self.gaussian_fit(t_refined, A, B, C)

                # Apply full width at half maximum (FWHM) and EIM intensity thresholds
                if self.fwhm_threshold(C, dt_i):
                    round_B = str(round(B, 2))
                    dt = float(round_B) if round_B else None

                # Convert extracted drift times to calibrated CCS values using CCSCalibrationRawXL module
                ccs = None
                if dt is not None:
                    ccs = self.cal_data.calibrated_ccs(monoisotopic_mz, dt)
                    ccs = round(ccs, 2)

                    # Check if sample is an internal standard
                    if sample_type == "IS":

                        # Create a key for looking up the reference CCS based on m/z (and gradient)
                        ref_ccs_key = (mz, gradient)

                        # If the key is in the reference dictionary, check against the reference CCS value
                        if ref_ccs_key in self.reference_is_dict:
                            reference_ccs = self.reference_is_dict[ref_ccs_key][
                                "Reference CCS (Å²)"
                            ]

                            # Calculate the percentage difference between observed and reference CCS
                            ccs_difference = abs(ccs - reference_ccs) / reference_ccs

                            # Flag the CCS if the difference is greater than 3%
                            if ccs_difference > 0.03:
                                ccs_flags[ccs] = "FLAG"
                else:
                    ccs = None

                # Calculate the EIC peak height (i.e., max intensity)
                peak_height = max(rt_i[start_idx : end_idx + 1])

                # Append the extracted data to output_rows list
                flag_status = ccs_flags.get(ccs, "")
                ccs_output = "" if ccs is None else f"{ccs} {flag_status}".strip()
                self.output_rows.append(
                    [
                        file_name,
                        compound_name,
                        adduct,
                        mz,
                        monoisotopic_mz,
                        sample_type,
                        row["CCS Calibrant"],
                        row["Gradient"],
                        row["Column Type"],
                        dt,
                        ccs_output,
                        rt_value,
                        peak_height,
                        peak_area,
                    ]
                )

                # Generate title for ATD figure
                fwhm = C * 2.355
                title_atd = f"m/z: {mz:.4f} \u00B1 {mz_tolerance} rt: {rt_value} → FWHM ~ {fwhm:.2f}"

                # Generate file name for EIM figure without the directory path
                fname_suffix = "IS_EIM" if sample_type == "IS" else "EIM"
                fname = "{}_{}_{}_{}.png".format(file_name, mz, rt_value, fname_suffix)

                # Replace spaces with underscores in the filename only
                fname = fname.replace(" ", "_")

                # Prepend the directory path to the filename
                """fname_atd = os.path.join(mobilogram_directory, fname)

                # Generate EIM figure
                self.atd(
                    t, t_refined, dt_i, fit_i, A, B, title_atd, fname_atd
                )"""  # Comment out this code if figures are not needed

                # Generate file name for combined figure
                fname_combined = (
                    f"{data_directory}/{compound_name}_{mz}_{adduct}_{file_name}.png"
                )

                # Generate combined figure
                self.combined_figure(
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
                    rt_value,
                    fwhm,
                    mu_values,
                    fname_combined,
                )

        # Export DataFrame containing extracted spectral features to Excel file (.xlsx)
        self.output_file = self.export_to_excel()

        # Export DataFrame containing filtered spectral features to Excel file (.xlsx)
        """self.filter_data()"""
