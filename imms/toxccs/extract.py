"""
    toxccs/extract.py
    Ryan Nguyen
    2/20/2024

    description:
    Module designed to handle the extraction and basic processing of full-scan LC-IM-MS/MS data. Provides a comprehensive set of tools for extracting raw data from .raw files, detecting peaks in the time and m/z dimensions, and converting drift times to calibrated collision cross-section (CCS) values. 
"""

import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, gaussian, convolve
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
import warnings
from dhrmasslynxapi.reader import MassLynxReader
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRawXL
from rdkit import Chem
from rdkit.Chem import Draw
import io
from PIL import Image

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

    def smiles_to_structure(self, smiles, img_size=(100, 100)):
        """
        RawProcessor.smiles_to_structure
        description:
                Converts a SMILES string to a chemical structure.
        parameters:
                smiles (str) -- SMILES string.
                img_size (tuple) -- image size. Default is (100, 100).
        """

        # Check if the string is valid
        if not isinstance(smiles, str) or smiles.strip() == "":
            return None
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            img = Draw.MolToImage(mol, size=img_size)
            return img
        else:
            return None

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
        RawProcessor.combined_figure
        description:
                Generates a figure containing the extracted chromatogram, mobilogram, and MS1 scan with the identified monoisotopic peak displayed.
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
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(6.4, 9.6))

        # Plot LC chromatogram figure for successfully fitted EICs
        ax2.plot(rt, rt_i, "b-", lw=1.5, label="Raw Data")
        ax2.plot(rt, smoothed_rt_i, "g--", lw=1.5, label="Smoothed Data")

        # Condition to plot if peaks were successfully fitted
        if len(peak_indices) > 0:
            ax2.plot(
                rt,
                self.multi_gaussian_fixed_mu(rt, mu_values, *popt_multi_gaussian),
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
        y_values = self.multi_gaussian_fixed_mu(
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
        ax2.set_yticklabels(
            [int(y) for y in y_tick_values], fontdict=tick_label_fontprops
        )
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
        chem_struct_img = self.smiles_to_structure(smiles, img_size=(350, 350))

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
        ax1.axvline(
            x=monoisotopic_mz, color="magenta", label="Monoisotopic Peak", lw=2.5
        )
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
                "SMILES",
                "Adduct",
                "Target m/z",
                "Observed m/z",
                "Sample Type",
                "CCS Calibrant",
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
                (tuple) -- extracted (m/z, intensity) values of the monoisotopic peak.
        """

        # Filter peaks within the m/z tolerance range
        potential_peaks = [
            peak
            for peak in identified_peaks
            if abs(peak[0] - target_mz) <= mz_tolerance
        ]

        # If no peaks are found within tolerance, return the target m/z
        if not potential_peaks:

            return (target_mz, 0), 0

        # Assume the most intense peak is the monoisotopic peak
        monoisotopic_peak = max(potential_peaks, key=lambda x: x[1])

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

        return monoisotopic_peak, isotopologue_type

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
        returns:
                (.xlsx) -- Excel spreadsheet containing the extracted data for each feature.
                (.png) -- Image displaying the chemical structure, MS1 scan, EIC, and EIM.
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

        # Generate a new folder called "Extracted Data" in the directory if not already present
        data_directory = "Extracted Data"
        if not os.path.exists(data_directory):
            os.makedirs(data_directory)

        # Iterate over each row in the target list DataFrame
        print(
            "\n...CCS calibration successful. Extracting data from the raw files using precursor target list {}...".format(
                self.target_list
            )
        )
        for i, row in df_input.iterrows():

            # Extract the path to .raw file, m/z, and sample type (i.e., sample or blank)
            file_name, mz, sample_type, gradient, adduct, compound_name, smiles = (
                row["File Name"],
                row["Target m/z"],
                row["Sample Type"],
                row["Gradient"],
                row["Adduct"],
                row["Compound Name"],
                row["SMILES"],
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
            (monoisotopic_mz, monoisotopic_intensity), _ = self.observed_mz(
                identified_peaks, mz, self.mz_tolerance
            )

            if monoisotopic_mz == 0:
                continue

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
                else:
                    ccs = None

                # Calculate the EIC peak height (i.e., max intensity)
                peak_height = max(rt_i[start_idx : end_idx + 1])

                # Append the extracted data to output_rows list
                self.output_rows.append(
                    [
                        file_name,
                        compound_name,
                        smiles,
                        adduct,
                        mz,
                        monoisotopic_mz,
                        sample_type,
                        row["CCS Calibrant"],
                        row["Gradient"],
                        row["Column Type"],
                        dt,
                        ccs,
                        rt_value,
                        peak_height,
                        peak_area,
                    ]
                )

                # Calculate FWHM for figure
                fwhm = C * 2.355

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
                    monoisotopic_intensity,
                    rt_value,
                    fwhm,
                    mu_values,
                    fname_combined,
                    mz_spectrum,
                    intensity_spectrum,
                    smiles,
                )

        # Export DataFrame containing extracted spectral features to Excel file (.xlsx)
        self.output_file = self.export_to_excel()
