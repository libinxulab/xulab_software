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
from multigauss import process_data
from filter import filter_data

# Set global font conditions for figures
params = {"font.size": 8,
         "font.family": "Arial",
          "font.weight": "bold"}
plt.rcParams.update(params)

class RawProcessor:
	def __init__(self, target_list):
		"""
		RawProcessor.__init__
		description:
			Initializes a new RawProcessor with a target list containing precursor ions and 
			paths to .raw data files.
		parameters:
			target_list (str) -- path to the Excel file (.xlsx) containing the target list.
		"""
		self.target_list = target_list
		self.output_rows = []

	def peak_fit(self, t, dt_i, p0="guess"):
		"""
		RawProcessor.peak_fit
		description:
			Fits a peak in the EIM to a Gaussian curve.
		parameters:
			t (array) -- drift times.
			dt_i (array) -- intensity  values corresponding to the drift times.
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
			t (array) -- raw EIM time points.
			t_refined (array) -- refined data points for plotting the fit.
			dt_i (array) -- original intensity values.
			fit_i (array) -- intensity values of the fitted curve.
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
			Exports the extracted data to an Excel (.xlsx) file.
		returns:
		(str) -- path to the output Excel (.xlsx) file.
		"""
		df_output = pd.DataFrame(self.output_rows, columns=["file_name", "mz", "sample_type", "ccs_calibrant", "gradient", "column_type", "monoisotopic_mz", "dt", "ccs", "rt", "peak_area"])
		output_file = self.target_list[:-5] + "_processed.xlsx"
		df_output.to_excel(output_file, index=False)
		print("\nProcessing successful. View results in {}".format(output_file))
		return output_file

	def filter_data(self):
		"""
		RawProcessor.filter_data
		description:
			Filters the extracted data to identify unique spectral features based on LC peak intensity. 
		returns:
		(str) -- The path to the filtered output Excel (.xlsx) file.
		"""
		# Read the input Excel file containing extracted data, i.e., the output of RawProcess.extract
		df = pd.read_excel(self.output_file)

		# Filter rows where sample_type = sample
		sample_rows = df[df["sample_type"] == "sample"]
		unique_features = []

		for _, sample_row in sample_rows.iterrows():
			sample_mz = sample_row["mz"]
			sample_rt = sample_row["rt"]
			sample_peak_area = sample_row["peak_area"]

			# Filter control rows with matching m/z and rt within +/- 0.1 min
			control_rows = df[(df["sample_type"] == "control") & (df["mz"] == sample_mz) & (abs(df["rt"] - sample_rt) <=0.1)]
			if len(control_rows) == 0:
				unique_features.append(sample_row)
			else:
				control_peak_area = control_rows["peak_area"].max()

				# Remove extracted features that have peak areas 1.5 times less than control features
				if sample_peak_area > 1.5 * control_peak_area:
					unique_features.append(sample_row)

		# Export filtered spectral features to pandas DataFrame
		df_output = pd.DataFrame(unique_features)
		filtered_file = self.output_file[:-15] + "_filtered.xlsx"

		# Export pandas DataFrame to Excel (.xlsx) file
		df_output.to_excel(filtered_file, index=False)
		print("\nFiltering sucessful. View results in {}\n".format(filtered_file))

	def monoisotopic_peak(self, mz_values, intensities, target_mz, tolerance=3):
		"""
		RawProcessor.monoisotopic_peak
		description:
			Identifies the monoisotopic peak within the target rt-selected MS1 spectrum.
			The current implementation only considers carbon isotopes around a M+3 window. 
			Method adapted from MetaboAnnotatoR (Ebbels et al., 2022)
		parameters:
			mz_values (array) -- extracted m/z values from the MS1 function.
			intensities (array) -- corresponding intensity values from the MS1 function.
			target_mz (float) -- target precursor m/z value to compare with extracted peaks.
		returns:
			(float) -- extracted m/z value of the monoisotopic peak. 
		"""
		# Identify the peaks from the extracted MS1 function that fall within the M+3 window
		potential_peaks = [(mz, intensity) for mz, intensity in zip(mz_values, intensities) if target_mz - tolerance <= mz <= target_mz]
		if not potential_peaks:
			return target_mz

		# Find the peak with the highest intensity within the tolerance window
		highest_intensity_peak = max(potential_peaks, key=lambda x: x[1])

		# If the target m/z has the highest intensity, it's assumed to be the monoisotopic peak
		# Otherwise, return the m/z of the highest intensity peak
		return target_mz if highest_intensity_peak[0] == target_mz else highest_intensity_peak[0]

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

		# Iterate over each row in the target list DataFrame
		print("\n...CCS Calibration successful. Extracting data from the raw files using precursor target list {}...".format(self.target_list))
		for i, row in df_input.iterrows():

			# Extract the path to .raw file, m/z, and sample type (i.e., sample or blank)
			file_name, mz, sample_type = row["file_name"], row["mz"], row["sample_type"]

			# Print the current m/z being queried to the terminal
			print("\n\traw file: {} m/z: {}".format(file_name, mz))

			# Extract m/z-selected ion chromatogram (EIC) from the .raw file using dhrmasslynxapi
			# Requires user-inputted MS1 function number and desired m/z tolerance
			rdr = MassLynxReader(file_name)
			rt, rt_i = rdr.get_chrom(self.ms1_function, float(mz), self.mz_tolerance)

			# Store extracted rt, rt_i data as arrays for further processing
			rt = np.array(rt)
			rt_i = np.array(rt_i)

			# Smooth and fit EIC with multigauss module
			peak_indices, rt_list, areas = process_data(rt, rt_i, file_name, mz, sample_type, generate_figures=True)

			# Iterate over each extracted LC ion chromatogram peak
			for rt_value, peak_area in zip(rt_list, areas):

				# Extract MS1 spectrum for isotopologue check
				mz_spectrum, intensity_spectrum = rdr.get_spectrum(ms1_function, float(rt_value) - 0.1, float(rt_value) + 0.1)

				# Identify the monoisotopic peak
				monoisotopic_mz = self.monoisotopic_peak(mz_spectrum, intensity_spectrum, float(mz))


				# Extract (m/z,rt)-selected ion mobilogram (EIM) for each identified LC peak using +/- 0.1 min window
				t, dt_i = rdr.get_filtered_chrom(self.mobility_function, float(mz), self.mz_tolerance, rt_min=float(rt_value)-0.1, rt_max=float(rt_value)+0.1)

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

				# Append the extracted data to output_rows list
				self.output_rows.append([file_name, mz, sample_type, row["ccs_calibrant"], row["gradient"], row["column_type"], monoisotopic_mz, dt, ccs, rt_value, peak_area])

				# Generate title for EIM figure
				fwhm = C * 2.355
				title_atd = "Extracted Mobilogram ({}) \nm/z: {:.4f} \u00B1 0.025 rt: {} → FWHM ~ {:.2f} ms".format(file_name, mz, rt_value, fwhm)

				# Generate file name for EIM figure
				fname_atd = "{}/{}_{}_{}_{}_EIM.png".format(mobilogram_directory, file_name, mz, rt_value, sample_type.capitalize()) # Change where mobilogram figures are saved
			
				# Generate EIM figure
				self.atd(t, t_refined, dt_i, fit_i, A, B, title_atd, fname_atd) # Comment out this code if figures are not needed
		
		# Export DataFrame containing extracted spectral features to Excel (.xlsx) file
		self.output_file = self.export_to_excel()

		# Export DataFrame containing filtered spectral features to Excel (.xlsx) file
		self.filter_data()
		 