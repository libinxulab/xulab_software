"""
toxccs.spectra.py
Ryan Nguyen
2/26/2024

description:
Module designed to handle the extraction and basic processing of fragmentation spectra from LC-IM-MS/MS data. Provides a comprehensive set of tools for extracting fragmentation spectra from .raw files, applying flters from the LC and mobility dimensions, and saving the data in convenient formats for downstream applications. 
"""

import pandas as pd
import sqlite3
import numpy as np
import math
import os
import re
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from scipy.signal import convolve, gaussian, find_peaks
from dhrmasslynxapi.reader import MassLynxReader
from rdkit import Chem
from rdkit.Chem import Draw
import io
from PIL import Image

# Set global font conditions for figures
params = {"font.family": "Arial", "font.weight": "bold"}
plt.rcParams.update(params)


class SpectralProcessor:
    def __init__(self, feature_list):
        """
        SpectralProcessor.__init__
        description:
                Initializes a new SpectralProcessor with a feature list containing precursor ions, retention times, mobility data, and paths to .raw data files.
        parameters:
                feature_list (str) -- path to the Excel file (.xlsx). containing the target list.
        """
        self.feature_list = feature_list

    def smiles_to_structure(self, smiles, img_size=(100, 100)):
        """
        SpectralProcessor.smiles_to_structure
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

    def export_to_excel(self, processed_spectrum, file_name, compound_name, adduct):
        """
        SpectralProcessor.export_to_excel
        description:
                Exports the processed spectrum data, including the original and normalied intensities, and top 5 most intense peaks to an Excel file.

        returns:
                (str) -- path to the Excel file (.xlsx).
        """
        # Prepare the DataFrame
        df = pd.DataFrame(
            {
                "m/z": processed_spectrum["mz"],
                "Intensity": processed_spectrum["intensity_raw"],
                "Normalized Intensity": processed_spectrum["normalized_intensity"],
                "Top 5 Peak m/z": None,
                "Top 5 Peak Intensity": None,
                "Top 5 Peak Normalized Intensity": None,
            }
        )

        # Prepare a DataFrame for the top 5 peaks
        top_5_peaks_df = pd.DataFrame(processed_spectrum["top_5_peak_info"])

        # Sort by m/z value
        top_5_peaks_df = top_5_peaks_df.sort_values(by="mz").reset_index(drop=True)

        # Fill in the top 5 peaks information
        for i in range(len(top_5_peaks_df)):
            df.loc[i, "Top 5 Peak m/z"] = top_5_peaks_df.loc[i, "mz"]
            df.loc[i, "Top 5 Peak Intensity"] = top_5_peaks_df.loc[i, "intensity_raw"]
            df.loc[i, "Top 5 Peak Normalized Intensity"] = top_5_peaks_df.loc[
                i, "intensity_normalized_to_top_5"
            ]

        # Define the output file path and save the DataFrame
        output_path = os.path.join(
            "Extracted Spectra", f"{compound_name}_{adduct}_{file_name}.xlsx"
        )
        df.to_excel(output_path, index=False)

        return output_path

    def peak_fit(self, t, dt_i, p0="guess"):
        """
        SpectralProcessor.peak_fit
        description:
                 Fits a peak in the EIM to a Gaussian curve.
        parameters:
                t (array) -- drift times.
                dt_i (array) -- intensity values corresponding to the drift times.
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
        SpectralProcessor.gaussian_fit
        description:
                Gaussian function for fitting EIMs.
        parameters:
                x (array) -- data points where Gaussian function is evaluated.
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
        SpectralProcessor.fwhm_thresholds
        description:
                Checks if the full width at half maximum (FWHM) and intensity of the EIM meet threshold values.
        parameters:
                C (float) -- standard deviation of the Gaussian curve.
                dt_i (array) -- intensity values.
                fwhm_thresholds (tuple) -- minimum and maximum acceptable FWHM.
                intensity_threshold (float) -- minimum acceptable intensity.
        returns:
                (bool) -- True if the FWHM and intensity meet the thresholds, False otherwise.
        """
        fwhm = C * 2.355
        max_intensity = max(dt_i)

        return (
            fwhm_thresholds[0] < fwhm < fwhm_thresholds[1]
            and max_intensity > intensity_threshold
        )

    def gaussian_smooth_pick(
        self, mz_array, intensity_array, window_len=1, std=0.1, prominence=0.001
    ):
        """
        SpectralProcessor.gaussian_smooth_pick
        description:
                Applies a Gaussian smoothing to the intensity array and identifies peaks.
                This function convolves the intensity array with a Gaussian window to smooth the data.
                It then identifies local maxima (i.e., peaks) in the smoothed data based on a defined prominence.
        parameters:
                mz_array (array-like): array of m/z values corresponding to the intensities.
                intensity_array (array-like): array of intensity values to be smoothed and from which peaks are identified.
                window_len (int) -- length of the gaussian window used for smoothing. Default is 1.
                std (float) -- standard deviation of the Gaussian window, controlling the degree of smoothing. Default is 0.1.
                prominence (float) -- minimum prominence of peaks to be identified. Default is 0.001.
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

    def extract_aif_rt(self, file_name, ms2_function, rt_start, rt_end):
        """
        SpectralProcessor.extract_aif_rt
        description:
                Extracts the raw retention time-selected AIF spectrum for a given spectral feature.
        parameters:
                file_name (str) -- path to the .raw file.
                ms2_function (int) -- MS2 function number.
                rt_start (float) -- retention time of the start of the identified EIC peak.
                rt_end (float) -- retention time of the end of the identified EIC peak.
        """
        # Initialize MassLynxReader object
        rdr = MassLynxReader(file_name)

        # Extract raw rt-selected AIF from the spectral feature using retention time indices
        m, i = rdr.get_spectrum(ms2_function, rt_start, rt_end)

        return np.array(m), np.array(i)

    def ms2_chromatogram(self, file_name, ms2_function, mz, mz_tolerance):
        """
        SpectralProcessor.ms2_chromatogram
        description:
                Extracts the MS2 level chromatogram for a given m/z and tolerance.
        parameters:
                file_name (str) -- path to the .raw file.
                ms2_function (int) -- MS2 function number.
                mz (float) -- m/z value used to extract the MS2 level chromatogram (presumably the observed monoisotopic peak).
                mz_tolerance (float) -- m/z tolerance for extraction.
        returns:
                (tuple): array of retention times and corresponding intensity values.
        """
        # Initialize MassLynxReader object
        rdr = MassLynxReader(file_name)

        # Extract the m/z-selected MS2 level chromatogram
        rt_array, rt_i_array = rdr.get_chrom(ms2_function, mz, mz_tolerance)

        return np.array(rt_array), np.array(rt_i_array)

    def extract_aif(
        self, file_name, mobility_function, rt_start, rt_end, dt_start, dt_end
    ):
        """
        SpectralProcessor.extract_aif
        description:
                Extracts the raw retention AND drift time-selected selected AIF spectrum for a given spectral feature.
        parameters:
                file_name (str) -- path to the .raw file.
                mobility_function (int) -- mobility function number.
                rt_start (float) -- retention time of the start of the identified EIC peak.
                rt_end (float) -- retention time of the end of the identified EIC peak.
                dt_start (float) -- drift time of the start of the identified EIM peak.
                dt_end (float) -- drift time of th end of the identified EIM peak.
        returns:
                (tuple) -- arrays of m/z values and corresponding intensities.
        """
        # Initialize MassLynxReader object
        rdr = MassLynxReader(file_name)

        # Extract raw (rt,dt)-selected AIF from the spectral feature using retention and drift time peak indices
        m, i = rdr.get_spectrum(
            mobility_function, rt_start, rt_end, dt_min=dt_start, dt_max=dt_end
        )

        return np.array(m), np.array(i)

    def extract_eim(
        self, file_name, function_number, mz, mz_tolerance, rt_start, rt_end
    ):
        """
        SpectralProcessor.extract_eim
        description:
                Extracts the raw retention time-selected EIM for a given m/z value.
        parameters:
                file_name (str) -- path to the .raw file.
                function_number (int) -- function number (i.e., MS1 or MS2 level, depending on status of the peak being processed).
                mz (float) -- m/z value used to extract the mobilogram.
                rt_start (float) -- start of the retention time range.
                rt_end (float) -- end of the retention time range.
        returns:
                (tuple) -- arrays of drift times and corresponding intensity vaules.
        """
        # Initialize MassLynxReader object
        rdr = MassLynxReader(file_name)

        # Extract m/z + retention time-selected mobilogram
        dt, intensity = rdr.get_filtered_chrom(
            function_number, float(mz), mz_tolerance, rt_start, rt_end
        )

        return np.array(dt), np.array(intensity)

    def process_extracted_spectrum(
        self, mz_array, intensity_array, monoisotopic_mz, mz_tolerance=0.025
    ):
        """
        SpectralProcessor.process_extracted_spectrum
        description:
                Processes the extracted AIF spectrum to normalize intensities and identify the top 5 most unique (i.e., non-isotopologue) peaks.
        parameters:
                mz_array (numpy.ndarray) -- array of m/z values.
                intensity_array (numpy.ndarray) -- array of intensity values.
                monoisotopic_mz (float) -- observed monoisotopic peak from MS1 scan used to remove non-relevant peaks.
                mz_tolerance (float) -- m/z tolerance around the monoisotopic peak for filtering. Default is 0.025 Da.
        returns:
                dict: a dictionary comtaining the raw intensities, normalized intensities, and the top 5 most intense peaks with their intensities normalized to the most intense peak among the top 5.

        """
        # Normalize the intensities to the most intense peak
        normalized_intensity = intensity_array / np.max(intensity_array) * 100

        # Remove peaks that are above the monoisotopic peak + 0.025 Da
        below_monoisotopic_indices = np.where(
            mz_array < monoisotopic_mz + mz_tolerance
        )[0]

        # Apply filter to mz, intensity, and normalized intensity arrays
        filtered_mz = mz_array[below_monoisotopic_indices]
        filtered_intensity = intensity_array[below_monoisotopic_indices]
        filtered_normalized_intensity = normalized_intensity[below_monoisotopic_indices]

        # Sort peaks by descending intensity
        sorted_indices = np.argsort(filtered_intensity)[::-1]
        sorted_mz = filtered_mz[sorted_indices]
        sorted_intensity = filtered_intensity[sorted_indices]
        sorted_normalized_intensity = filtered_normalized_intensity[sorted_indices]

        # Initialize the list to store indices of the top 5 peaks
        top_5_indices = []
        for i in range(len(sorted_mz)):
            if len(top_5_indices) >= 5:

                # Stop when we have collected the top 5 unique peaks
                break
            current_mz = sorted_mz[i]

            # Ensure the current peak is not within M+3 of any peak already in the top 5
            if not any(abs(current_mz - sorted_mz[j]) <= 3 for j in top_5_indices):
                top_5_indices.append(i)

        # Extract the final top 5 unique peaks using the indices
        final_top_5_mz = sorted_mz[top_5_indices]
        final_top_5_intensity_raw = sorted_intensity[top_5_indices]
        final_top_5_normalized_intensity = sorted_normalized_intensity[top_5_indices]

        # Construct the top 5 peak information
        top_5_peak_info = [
            {"mz": mz, "intensity_raw": raw, "intensity_normalized_to_top_5": norm}
            for mz, raw, norm in zip(
                final_top_5_mz,
                final_top_5_intensity_raw,
                final_top_5_normalized_intensity,
            )
        ]

        return {
            "mz": mz_array,
            "intensity_raw": intensity_array,
            "normalized_intensity": normalized_intensity,
            "top_5_peak_info": top_5_peak_info,
        }

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
        SpectralProcessor.process_chromatogram
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
                -- peak_ranges -- start and end indices for each identified peak.
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

    def combined_figure(
        self,
        mz_spectrum,
        intensity_spectrum,
        top_5_peak_info,
        observed_mz,
        compound_name,
        adduct,
        smiles,
        fname_combined,
        mz,
        rt_ms1,
        rt_i_ms1,
        rt_ms2,
        rt_i_ms2,
        eim_data,
    ):
        """
        SpectralProcessor.combined_figure
        description:
                Generates a figure containing the extracted AIF spectrum and a separate spectrum showing the top 5 most intense peaks.
        parameters:
                mz_spectrum (numpy.ndarray) -- array of m/z values.
                intensity_spectrum (numpy.ndarray) -- array of intensity values.
                top_5_peak_info (dict) -- dictionary containing m/z, raw intensity, and normalized intensity values of the top 5 most intense peaks from the AIF spectrum.
                observed_mz (float) -- observed monoisotopic m/z value (i.e., the parent ion).
                compound_name (str) -- name of the compound being analyzed.
                adduct (str) -- adduct type of the compound being analyzed.
                smiles (str) -- SMILES string of the compound being analyzed.
                fname_combined (str) -- file name for saving the plot.
                mz (float) -- target m/z value used to extract the raw data.
        """
        # Plot setup
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, figsize=(6.4, 9.6))

        # Plot full MS2 spectrum
        ax1.plot(mz_spectrum, intensity_spectrum, "k-", lw=1, label="MS2 Data")

        # Annotate the top 5 most intense peaks
        intensity_range = max(intensity_spectrum) - min(intensity_spectrum)
        normalized_offset_percentage = 0.05
        normalized_offset = intensity_range * normalized_offset_percentage
        for i, peak in enumerate(top_5_peak_info):

            # Draw vertical lines matching peak intensity
            if i == 0:
                ax1.vlines(
                    x=peak["mz"],
                    ymin=0,
                    ymax=peak["intensity_raw"],
                    colors="magenta",
                    lw=2,
                    label="Top 5 Peaks",
                )
            else:
                ax1.vlines(
                    x=peak["mz"],
                    ymin=0,
                    ymax=peak["intensity_raw"],
                    colors="magenta",
                    lw=2,
                )
            adjusted_x_position = peak["mz"]
            adjusted_y_position = peak["intensity_raw"] + normalized_offset
            ax1.text(
                adjusted_x_position,
                adjusted_y_position,
                f"{peak['mz']:.4f}",
                ha="center",
                va="top",
                fontsize=10,
                color="black",
                fontweight="bold",
                fontname="Arial",
            )
        ax1.legend(
            loc="best", fontsize=10, frameon=True, edgecolor="black", facecolor="white"
        )
        # Format MS2 spectrum
        ax1.set_xlim(50, observed_mz + 10)
        max_intensity_y_limit = max(intensity_spectrum) + 0.1 * max(intensity_spectrum)
        y_tick_values = np.linspace(0, max_intensity_y_limit, num=10, endpoint=True)
        ax1.set_yticks(y_tick_values)
        tick_label_fontprops = {"weight": "bold", "family": "Arial", "size": 10}
        ax1.set_yticklabels(
            [int(y) for y in y_tick_values], fontdict=tick_label_fontprops
        )
        ax1.set_ylim(0, max_intensity_y_limit)
        ax1.set_xlabel("m/z", fontweight="bold", fontname="Arial", fontsize=10)
        ax1.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
        ax1.set_title(
            f"{compound_name}\n{mz:.4f} {adduct}\nMS2 Spectrum",
            fontsize=12,
            fontweight="bold",
            fontname="Arial",
        )
        legend1 = ax1.legend(
            loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white"
        )
        for text in legend1.get_texts():
            text.set_fontname("Arial")
        ax1.spines["top"].set_visible(False)
        ax1.spines["right"].set_visible(False)
        ax1.spines["left"].set_linewidth(1.5)
        ax1.spines["bottom"].set_linewidth(1.5)
        ax1.tick_params(axis="both", which="major", labelsize=10, width=1)
        for label in ax1.get_xticklabels():
            label.set_fontname("Arial")
            label.set_fontweight("bold")
            label.set_fontsize(10)

        # Plot MS1 and MS2 level chromatograms
        ax2.plot(rt_ms1, rt_i_ms1, "b-", lw=1.75, label="MS1 Chromatogram")
        ax2.plot(rt_ms2, rt_i_ms2, "r-", lw=1.75, label="MS2 Chromatogram")
        ax2.set_xlabel(
            "Retention Time [min]", fontweight="bold", fontname="Arial", fontsize=10
        )
        ax2.set_xlim(0, 1.3)
        ax2.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
        ax2.set_title(
            "Precursor Ion Extracted Chromatograms",
            fontsize=12,
            fontweight="bold",
            fontname="Arial",
        )
        legend2 = ax2.legend(
            loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white"
        )
        for text in legend2.get_texts():
            text.set_fontname("Arial")
        max_intensity_y_limit = max(rt_i_ms1) + 0.1 * max(rt_i_ms1)
        y_tick_values = np.linspace(0, max_intensity_y_limit, num=10, endpoint=True)
        ax2.set_yticks(y_tick_values)
        ax2.set_yticklabels(
            [int(y) for y in y_tick_values], fontdict=tick_label_fontprops
        )
        ax2.set_ylim(0, max_intensity_y_limit)
        ax2.spines["top"].set_visible(False)
        ax2.spines["right"].set_visible(False)
        ax2.spines["left"].set_linewidth(1.5)
        ax2.spines["bottom"].set_linewidth(1.5)
        ax2.tick_params(axis="both", which="major", labelsize=10, width=1)
        for label in ax2.get_xticklabels():
            label.set_fontname("Arial")
            label.set_fontweight("bold")
            label.set_fontsize(10)

        # Plot EIMs
        latest_drift_time_apex = 0
        colors = ["red", "green", "blue", "orange", "purple"]
        for i, ((dt, intensity), peak_info) in enumerate(
            zip(eim_data, top_5_peak_info)
        ):
            ax3.plot(
                dt,
                intensity,
                color=colors[i],
                lw=1.75,
                label=f"{peak_info['mz']:.4f} m/z",
            )
            apex_index = np.argmax(intensity)
            apex_drift_time = dt[apex_index]
            latest_drift_time_apex = max(latest_drift_time_apex, apex_drift_time)
        ax3.set_xlabel(
            "Drift Time [ms]", fontsize=10, fontweight="bold", fontname="Arial"
        )
        max_intensity_all_mobilograms = 0
        for _, intensity in eim_data:
            max_intensity = max(intensity)
            max_intensity_all_mobilograms = max(
                max_intensity_all_mobilograms, max_intensity
            )
        max_intensity_y_limit = (
            max_intensity_all_mobilograms + 0.1 * max_intensity_all_mobilograms
        )
        y_tick_values = np.linspace(0, max_intensity_y_limit, num=10, endpoint=True)
        ax3.set_yticks(y_tick_values)
        ax3.set_yticklabels(
            [int(y) for y in y_tick_values], fontdict=tick_label_fontprops
        )
        ax3.set_ylim(0, max_intensity_y_limit)
        ax3.set_ylabel("Intensity", fontsize=10, fontweight="bold", fontname="Arial")
        ax3.set_title(
            "Extracted Mobilograms", fontsize=12, fontweight="bold", fontname="Arial"
        )
        ax3.set_xlim(0, latest_drift_time_apex + 2)

        legend3 = ax3.legend(
            loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white"
        )
        for text in legend3.get_texts():
            text.set_fontname("Arial")
        ax3.tick_params(axis="both", which="major", labelsize=10, width=1)
        for label in ax3.get_xticklabels():
            label.set_fontname("Arial")
            label.set_fontweight("bold")
            label.set_fontsize(10)
        ax3.spines["top"].set_visible(False)
        ax3.spines["right"].set_visible(False)
        ax3.spines["left"].set_linewidth(1.5)
        ax3.spines["bottom"].set_linewidth(1.5)
        plt.tight_layout()
        plt.savefig(fname_combined, dpi=300, bbox_inches="tight")
        plt.close()

    def extract_spectra(
        self, ms1_function, ms2_function, mobility_function, mz_tolerance
    ):
        """
        SpectralProcessor.extract_spectra
        description:
                Main function for extracting fragmentation spectra from .raw data files. The observed monoisotpic mass of the feature is first used to extract the MS1 chromatogram bounds, which are then applied to extract drift time bounds. These bounds are then simultaneously applied to extract the (rt,dt)-selected fragmentation spectrum from the MS2 level.
        parameters:
                ms1_function (int) -- MS1 function number.
                ms2_function (int) -- MS2 function number.
                mobility_function (int) -- mobility function number.
        """
        # Set the input parameters
        self.ms1_function = ms1_function
        self.ms2_function = ms2_function
        self.mobility_function = mobility_function
        self.mz_tolerance = mz_tolerance

        # Read feature list into a pandas DataFrame
        df_input = pd.read_excel(self.feature_list)
        df_input.reset_index(drop=True, inplace=True)

        # Generate a new folder called "Extracted Spectra" in the directory if not already present
        data_directory = "Extracted Spectra"
        if not os.path.exists(data_directory):
            os.makedirs(data_directory)

        # Generate a new folder called "Extracted Spectra Figures" in the directory if not already present
        figures_directory = "Extracted Spectra Figures"
        if not os.path.exists(figures_directory):
            os.makedirs(figures_directory)

        # Iterate over each row in the feature list DataFrame
        print(
            "\n...Extracting fragmentation spectra for spectral features in {}...".format(
                self.feature_list
            )
        )

        for i, row in df_input.iterrows():

            # Extract the path to the .raw file, compound name, adduct type, observed monoisotopic m/z, retention time, drift time, sample type, and target m/z
            (
                file_name,
                compound_name,
                adduct,
                monoisotopic_mz,
                rt,
                dt,
                sample_type,
                mz,
                smiles,
            ) = (
                row["File Name"],
                row["Compound Name"],
                row["Adduct"],
                row["Observed m/z"],
                row["Observed Retention Time (min)"],
                row["Observed Drift Time (ms)"],
                row["Sample Type"],
                row["Target m/z"],
                row["SMILES"],
            )

            # Print the current feature being processed to the terminal
            print("\n\traw file: {} compound: {}".format(file_name, compound_name))

            # Extract m/z-selected ion chromatogram (EIC) from the .raw file using observed m/z value
            rdr = MassLynxReader(file_name)
            rt_array, rt_i_array = rdr.get_chrom(
                self.ms1_function, float(monoisotopic_mz), self.mz_tolerance
            )

            # Store extracted rt, rt_i data as arrays for further processing
            rt_array, rt_i_array = np.array(rt_array), np.array(rt_i_array)

            # Smooth and fit EIC
            (
                peak_indices,
                rt_list,
                areas,
                peak_ranges,
                popt_multi_gaussian,
                mu_values,
            ) = self.process_chromatogram(
                rt_array,
                rt_i_array,
                file_name,
                compound_name,
                adduct,
                mz,
                sample_type,
                generate_images=False,
            )

            # find the LC peak that is closest to the feature list value
            closest_peak_idx = np.argmin(np.abs(rt_array[peak_indices] - rt))

            # Find the LC peak range (i.e., the start and end indices)
            start_idx, end_idx = peak_ranges[closest_peak_idx]

            # Use the peak indices to define the EIC window
            rt_start, rt_end = rt_array[start_idx], rt_array[end_idx]

            # Handle cases where peak indices are determined to be equal due to close spacing of raw data
            # Default to use a fixed 0.3 min window (i.e., +/- 0.15 min around the apex)
            if rt_start == rt_end:
                fixed_bound = 0.15

                # Adjust rt_start and rt_end for these cases
                rt_start = max(rt[0], rt[start_idx] - fixed_bound)
                rt_end = min(rt[-1], rt[end_idx] + fixed_bound)

            # Use the EIC peak indices to extract the (m/z,rt)-selected ion mobilogram (EIM)
            # Need to use MS1 level mobility function number here
            t, dt_i = rdr.get_filtered_chrom(
                0,
                float(monoisotopic_mz),
                self.mz_tolerance,
                rt_min=rt_start,
                rt_max=rt_end,
            )

            # Smooth and fit EIM to Gaussian function
            A, B, C = self.peak_fit(t, dt_i)
            t_refined = np.arange(min(t), max(t), 0.01)
            fit_i = self.gaussian_fit(t_refined, A, B, C)
            fitted_dt = float(round(B, 2))

            # Ensure that the extracted drift time is equal to the drift time of the extracted spectral feature
            # Because this value is extracted using the same functions, it should be identical
            # Using a threshold just in case
            if abs(fitted_dt - dt) <= 0.2:

                # Calculate full width at 1% maximum (FW1M) bounds of EIM
                fw1m = 2 * np.sqrt(-2 * C**2 * np.log(0.01))

                # Find the EIM peak range (i.e., the start and end drift times corresponding to the FW1M)
                dt_start = B - (fw1m / 2)
                dt_end = B + (fw1m / 2)

                # Ensure dt_start and dt_end are within the bounds of the original drift time array
                dt_start = max(dt_start, min(t))
                dt_end = min(dt_end, max(t))

            """# Extract rt-selected AIF spectrum for the feature using identified EIC bounds
            extracted_mz, extracted_intensity = self.extract_aif_rt(
                file_name, self.ms2_function, rt_start, rt_end
            )"""

            # Extract (rt,dt)-selected AIF spectrum for the feature using identified EIC and EIM bounds
            extracted_mz, extracted_intensity = self.extract_aif(
                file_name,
                self.mobility_function,
                rt_start,
                rt_end,
                dt_start,
                dt_end,
            )

            # Store and process the extracted AIF spectrum to get normalized intensities and top 5 most intense peaks
            processed_spectrum = self.process_extracted_spectrum(
                extracted_mz, extracted_intensity, monoisotopic_mz, self.mz_tolerance
            )

            # Export the processed spectrum to an Excel file for each compound
            output_excel_path = self.export_to_excel(
                processed_spectrum, file_name, compound_name, adduct
            )

            # Generate file name for combined figure
            fname_combined = (
                f"{figures_directory}/{compound_name}_{mz}_{adduct}_{file_name}.png"
            )

            # Fetch MS2 level chromatogram for combined figure
            rt_ms2, rt_i_ms2 = self.ms2_chromatogram(
                file_name, self.ms2_function, monoisotopic_mz, self.mz_tolerance
            )

            # Fetch EIMs for each of the top 5 most intense peaks in the extracted AIF spectrum
            eim_data = []
            for peak in processed_spectrum["top_5_peak_info"]:
                mz_peak = peak["mz"]
                if abs(mz_peak - monoisotopic_mz) <= self.mz_tolerance:

                    # If precursor peak is present, extract from MS1 level
                    dt, intensity = self.extract_eim(
                        file_name,
                        self.ms1_function,
                        mz_peak,
                        self.mz_tolerance,
                        rt_start,
                        rt_end,
                    )
                else:

                    # Extract from MS2 level
                    dt, intensity = self.extract_eim(
                        file_name,
                        self.ms2_function,
                        mz_peak,
                        self.mz_tolerance,
                        rt_start,
                        rt_end,
                    )
                eim_data.append((dt, intensity))

            # Generate combined figure
            self.combined_figure(
                mz_spectrum=extracted_mz,
                intensity_spectrum=extracted_intensity,
                top_5_peak_info=processed_spectrum["top_5_peak_info"],
                observed_mz=monoisotopic_mz,
                compound_name=compound_name,
                adduct=adduct,
                smiles=smiles,
                fname_combined=fname_combined,
                mz=mz,
                rt_ms1=rt_array,
                rt_i_ms1=rt_i_array,
                rt_ms2=rt_ms2,
                rt_i_ms2=rt_i_ms2,
                eim_data=eim_data,
            )
