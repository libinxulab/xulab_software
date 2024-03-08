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

# TODO
# Implement plotting diagram that displays the EIC and EIM + indices, full extracted spectrum, and top 5 peaks

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

    def extract_aif(
        self, file_name, mobility_function, rt_start, rt_end, dt_start, dt_end
    ):
        """
        SpectralProcessor.extract_aif
        description:
                Extracts the raw AIF spectrum for a given spectral feature.
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

    def process_extracted_spectrum(self, mz_array, intensity_array):
        """
        SpectralProcessor.process_extracted_spectrum
        decription:
                Processes the extracted AIF spectrum to normalize intensities and identify the top 5 most intense peaks.
        parameters:
                mz_array (numpy.ndarray): array of m/z values.
                intensity_array (numpy.ndarray): array of intensity values.
        returns:
                dict: a dictionary containing the raw intensities, normalized intensities, and the top 5 most intense peaks with their intensities normalized to the most intense peak among the top 5.
        """
        # Normalize the intensities to the most intense peak
        max_intensity_entire_spectrum = np.max(intensity_array)
        normalized_intensity_entire_spectrum = (
            intensity_array / max_intensity_entire_spectrum
        ) * 100

        # Identify the top 5 most intense peaks
        peak_indices = np.argsort(intensity_array)[::-1][:5]
        top_5_mz = mz_array[peak_indices]
        top_5_intensity_raw = intensity_array[peak_indices]

        # Normalize the top 5 most intense peaks to the most intense peak among them
        max_top_5_intensity = np.max(top_5_intensity_raw)
        normalized_top_5_intensity = top_5_intensity_raw / max_top_5_intensity * 100

        # Construct the top 5 peak information
        top_5_peak_info = [
            {"mz": mz, "intensity_raw": raw, "intensity_normalized_to_top_5": norm}
            for mz, raw, norm in zip(
                top_5_mz, top_5_intensity_raw, normalized_top_5_intensity
            )
        ]

        return {
            "mz": mz_array,
            "intensity_raw": intensity_array,
            "normalized_intensity": normalized_intensity_entire_spectrum,
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
            ) = (
                row["File Name"],
                row["Compound Name"],
                row["Adduct"],
                row["Observed m/z"],
                row["Observed Retention Time (min)"],
                row["Observed Drift Time (ms)"],
                row["Sample Type"],
                row["Target m/z"],
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
            t, dt_i = rdr.get_filtered_chrom(
                self.mobility_function,
                float(monoisotopic_mz),
                self.mz_tolerance,
                rt_min=rt_start,
                rt_max=rt_end,
            )

            # Parameter for smoothing Gaussian function
            t_refined = np.arange(min(t), max(t), float(0.001))

            # Initialzie and fit Gaussian function
            A, B, C = self.peak_fit(t, dt_i)
            fit_i = self.gaussian_fit(t_refined, A, B, C)
            fitted_dt = float(round(B, 2))

            # Ensure that the extracted drift time is equal to the drift time of the extracted spectral feature
            # Because this value is extracted using the same functions, it should be identical
            # Using a threshold just in case
            if abs(fitted_dt - dt) <= 0.1:

                # Calculate full width at 1% maximum (FW1M) bounds of EIM
                fw1m = 2 * np.sqrt(-2 * C**2 * np.log(0.01))

                # Find the EIM peak range (i.e., the start and end drift times corresponding to the FW1M)
                dt_start = B - (fw1m / 2)
                dt_end = B + (fw1m / 2)

                # Ensure dt_start and dt_end are within the bounds of the original drift time array
                dt_start = max(dt_start, min(t))
                dt_end = min(dt_end, max(t))

            # Extract (rt,dt)-selected AIF spectrum for the feature using identified EIC and EIM bounds
            extracted_mz, extracted_intensity = self.extract_aif(
                file_name, self.mobility_function, rt_start, rt_end, dt_start, dt_end
            )

            # Store and process the extracted AIF spectrum to get normalized intensities and top 5 most intense peaks
            processed_spectrum = self.process_extracted_spectrum(
                extracted_mz, extracted_intensity
            )

            # Export the processed spectrum to an Excel file for each compound
            output_excel_path = self.export_to_excel(
                processed_spectrum, file_name, compound_name, adduct
            )
