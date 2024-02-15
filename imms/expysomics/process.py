"""
	expysomics/process.py
	Ryan Nguyen
	2/15/2024

	description:
		Module designed to handle the extraction and basic processing of raw LC-IM-MS data. 
		Provides a comprehensive set of tools for extracting data from .raw files,
		peak detection in both the time and m/z dimensions, Gaussian fittin, and automatic 
		conversion of drift times to calibrated collision cross-section (CCS) values.
		This module may be used on its own to extract data or in combination with expysomics/query.py for 
		downstream analysis.
"""
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, gaussian, convolve
import math
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from dhrmasslynxapi.reader import MassLynxReader
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRawXL
from multigauss import process_chromatogram

# Set global font conditions for figures
params = {"font.size": 8,
		 "font.family": "Arial",
		  "font.weight": "bold"}
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

		return A * np.exp(-(x - B) ** 2 / (2 * C ** 2))

	def fwhm_threshold(self, C, dt_i, fwhm_thresholds=(0.05, 2.5), intensity_threshold=500):
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

		return fwhm_thresholds[0] < fwhm < fwhm_thresholds[1] and max_intensity > intensity_threshold

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
		fig, ax = plt.subplots(figsize=(6.4, 4.8)) 
		ax.plot(t, dt_i, "o--b", lw=1.5, ms=2, label="Raw Data") # Plot raw EIM
		ax.plot(t_refined, fit_i, "r-", lw=1.5, label="Gaussian Fit") # Plot fitted data
		legend = ax.legend(loc="best", frameon=True, fontsize=10, edgecolor="black", facecolor="white")
		for text in legend.get_texts():
			text.set_fontname("Arial")
		ax.text(B + 0.1, 0.95 * max(fit_i), "{:.2f}".format(B), c="k", fontsize=10, fontweight="bold", fontname="Arial") # Annotate peak of fitted data with dt
		ax.set_title(title_atd, fontsize=12, fontweight="bold", fontname="Arial")
		ax.set_xlabel("Drift Time [ms]", fontsize=12, fontweight="bold", fontname="Arial")
		ax.set_ylabel("Intensity", fontsize=12, fontweight="bold", fontname="Arial")
		max_dt_i_y_limit = max(dt_i) + 0.1 * max(dt_i)
		plt.ylim(0, max_dt_i_y_limit) # Set y axis range
		plt.xlim(0, B + 2) # Set x axis range
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

	def export_to_excel(self):
		"""
		RawProcessor.export_to_excel
		description:
			Exports the extracted data to an Excel file (.xlsx). 
		returns:
		(str) -- path to the output Excel file (.xlsx).
		"""
		df_output = pd.DataFrame(self.output_rows, columns=["File Name", "Sample Type", "Gradient", "Column Type", "Target m/z", "Observed m/z", "CCS Calibrant",  "Observed Drift Time (ms)", "Observed CCS (Å²)", "Observed Retention Time (min)", "EIC Peak Intensity", "EIC Peak Area"])
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
			self.reference_is_dict[key] = {"Reference Retention Time (min)": row["Reference Retention Time (min)"], "Reference CCS (Å²)": row["Reference CCS (Å²)"]}

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
		df["Observed Retention Time (min)"] = df["Observed Retention Time (min)"].apply(lambda x: float(x.split()[0]) if isinstance(x, str) else x)

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
				matching_controls = matching_mz_controls[abs(matching_mz_controls["Observed Retention Time (min)"] - sample_rt) <= 0.1]

			# Handle case where there are no control features that match by m/z and retention time, or if the sample peak area is greater than the max control feature peak 
			if matching_controls.empty or (sample_peak_area > 1.5 * matching_controls["EIC Peak Area"].max()):

				# Assign a peak area ratio if matching control features are found
				if not matching_controls.empty:
					sample_row["Control EIC Peak Area Ratio"] = sample_peak_area / matching_controls["EIC Peak Area"].max()
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

	def gaussian_smooth_pick(self, mz_array, intensity_array, window_len=1, std=7, prominence=0.1):
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
			std (float) -- standard deviation of the Gaussian window, controlling the degree of smoothing. Default is 7.
			prominence (float) -- minimum prominence of peaks to be identified. Default is 0.1
		returns:
			(tuple) -- a tuple containing identified peaks and smoothed intensity array.
		"""
		window = gaussian(window_len, std=std)
		smoothed_intensity = convolve(intensity_array, window/window.sum(), mode="same")
		peaks, _ = find_peaks(smoothed_intensity, prominence=prominence)
		identified_peaks = [(mz_array[p], smoothed_intensity[p]) for p in peaks]

		return identified_peaks, smoothed_intensity

	def monoisotopic_peak(self, identified_peaks, target_mz, tolerance=0.025):
		"""
		RawProcessor.monoisotopic_peak
		description:
			Identifies the monoisotopic peak within the target rt-selected MS1 spectrum.
			The current implementation only considers the M+0 peak because all QACs carry 
			a permanent positive charge.
		parameters:
			identified_peaks (tuple) -- extracted and "centroided" MS1 peaks 
			target_mz (float) -- target precursor m/z value to compare with extracted peaks.
			tolerance (float) -- tolerance around the target m/z vakue to consider.
			(float) -- extracted m/z value of the monoisotopic peak. 
		"""
		# Identify the peaks from the extracted MS1 function that fall within the M+3 window
		potential_peaks = [(mz, intensity) for mz, intensity in identified_peaks if target_mz - tolerance <= mz <= target_mz + tolerance]
		if not potential_peaks:
			return target_mz

		# Find the peak with the highest intensity within the tolerance window
		highest_intensity_peak = max(potential_peaks, key=lambda x: x[1])

		return highest_intensity_peak[0]

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
		mobilogram_directory = "Extracted Mobilograms"
		if not os.path.exists(mobilogram_directory):
			os.makedirs(mobilogram_directory)

		# Initialize dictionary for internal standard CCS flags
		ccs_flags = {}

		# Iterate over each row in the target list DataFrame
		print("\n...CCS Calibration successful. Extracting data from the raw files using precursor target list {}...".format(self.target_list))
		for i, row in df_input.iterrows():

			# Extract the path to .raw file, m/z, and sample type (i.e., sample or blank)
			file_name, mz, sample_type, gradient = row["File Name"], row["Target m/z"], row["Sample Type"], row["Gradient"]

			# Print the current m/z being queried to the terminal
			print("\n\traw file: {} m/z: {}".format(file_name, mz))

			# Extract m/z-selected ion chromatogram (EIC) from the .raw file using theoretical m/z value
			# Requires user-inputted MS1 function number and desired m/z tolerance
			rdr = MassLynxReader(file_name)
			rt, rt_i = rdr.get_chrom(self.ms1_function, float(mz), self.mz_tolerance)

			# Store extracted rt, rt_i data as arrays for further processing
			rt = np.array(rt)
			rt_i = np.array(rt_i)

			# Smooth and fit EIC with multigauss module
			peak_indices, rt_list, areas, peak_ranges = process_chromatogram(rt, rt_i, file_name, mz, sample_type, generate_images=True)

			# Iterate over each extracted LC ion chromatogram peak
			for rt_value, peak_area, (start_idx, end_idx) in zip(rt_list, areas, peak_ranges):
				if sample_type =="IS":

					# Create a key for looking up the reference retention time based on m/z and gradient type
					ref_rt_key = (mz, gradient)

					# If the key is in the reference dictionary, check against the reference retention time
					if ref_rt_key in self.reference_is_dict:
						reference_rt = self.reference_is_dict[ref_rt_key]["Reference Retention Time (min)"]

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

				# Extract MS1 spectrum for isotopologue check
				mz_spectrum, intensity_spectrum = rdr.get_spectrum(ms1_function, rt_start, rt_end)

				# Smooth and fit the extracted MS1 spectrum using Gaussian convolution 
				# This process is analogous to "centroiding" the profile data as outlined in MSnbase
				# The default parmaeters for the smoothing and picking (i.e., window_len=1, std=7, prominence=0.01) select the most intense peak for each ion distribution
				# More aggressive smoothing will lead to slightly different m/z peaks being picked
				identified_peaks, smoothed_intensity = self.gaussian_smooth_pick(mz_spectrum, intensity_spectrum)

				# Identify the monoisotopic peak in the MS1 centroided data
				monoisotopic_mz = self.monoisotopic_peak(identified_peaks, float(mz))

				# Extract (m/z,rt)-selected ion mobilogram (EIM) for each identified LC peak 
				t, dt_i = rdr.get_filtered_chrom(self.mobility_function, float(mz), self.mz_tolerance, rt_min=rt_start, rt_max=rt_end)
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
					ccs = self.cal_data.calibrated_ccs(mz, dt)
					ccs = round(ccs, 2)

					# Check if sample is an internal standard
					if sample_type == "IS":

						# Create a key for looking up the reference CCS based on m/z (and gradient)
						ref_ccs_key = (mz, gradient)

						# If the key is in the reference dictionary, check against the reference CCS value
						if ref_ccs_key in self.reference_is_dict:
							reference_ccs = self.reference_is_dict[ref_ccs_key]["Reference CCS (Å²)"]

							# Calculate the percentage difference between observed and reference CCS
							ccs_difference = abs(ccs - reference_ccs) / reference_ccs

							# Flag the CCS if the difference is greater than 3%
							if ccs_difference > 0.03:
								ccs_flags[ccs] = "FLAG"
				else:
					ccs = None

				# Calculate the EIC peak height (i.e., max intensity)
				peak_height = max(rt_i[start_idx:end_idx+1])

				# Append the extracted data to output_rows list
				flag_status = ccs_flags.get(ccs, "")
				ccs_output = "" if ccs is None else f"{ccs} {flag_status}".strip()
				self.output_rows.append([file_name, sample_type, row["Gradient"], row["Column Type"], mz, monoisotopic_mz, row["CCS Calibrant"], dt, ccs_output, rt_value, peak_height, peak_area])

				# Generate title for EIM figure
				fwhm = C * 2.355
				title_atd = "Extracted Mobilogram ({}) \nm/z: {:.4f} \u00B1 0.025 rt: {} → FWHM ~ {:.2f} ms".format(file_name, mz, rt_value, fwhm)

				# Generate file name for EIM figure without the directory path
				fname_suffix = "IS_EIM" if sample_type == "IS" else "EIM"
				fname = "{}_{}_{}_{}.png".format(file_name, mz, rt_value, fname_suffix) 

				# Replace spaces with underscores in the filename only
				fname = fname.replace(" ", "_")

				# Prepend the directory path to the filename
				fname_atd = os.path.join(mobilogram_directory, fname)
			
				# Generate EIM figure
				self.atd(t, t_refined, dt_i, fit_i, A, B, title_atd, fname_atd) # Comment out this code if figures are not needed
		
		# Export DataFrame containing extracted spectral features to Excel file (.xlsx) 
		self.output_file = self.export_to_excel()

		# Export DataFrame containing filtered spectral features to Excel file (.xlsx) 
		self.filter_data()
		 