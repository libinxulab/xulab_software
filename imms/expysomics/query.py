"""
	expysomics/query.py
	Ryan Nguyen

	description:
		Module designed for the annotation of spectral features.
		Spectral features identified in expysomics/process.py are matched to potential compounds
		based on various selection criteria (i.e., m/z, retention time, CCS, and/or fragmentation spectra).
		Potential matches are ranked according to various scores (i.e., cosine similarity, mass error, CCS error).
		This module is currently designed to query SQLite CCS and fragmentation reference databases and 
		accepts the output of expysomics/process.py.
"""
import pandas as pd 
import sqlite3
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from scipy.signal import convolve, gaussian, find_peaks
from multigauss import process_chromatogram
from dhrmasslynxapi.reader import MassLynxReader

# Set global font conditions for figures
params = {"font.family": "Arial",
		  "font.weight": "bold"}
plt.rcParams.update(params)

class FeatureAnnotate:
	def __init__(self, feature_list, reference_db, spectral_db=None):
		"""
		FeatureAnnotate.__init__
		description:
			Initializes a new FeatureAnnotate with a spectral feature list to be annotated.
		parameters:
			feature_list (str) -- path to the Excel file (.xlsx) containing the feature list.
			reference_db (str) -- path to the reference database file (.db)
			spectral_db (str) -- path to the reference spectral database file (.db)
		"""
		self.feature_list = feature_list
		self.feature_df = None
		self.reference_db = reference_db 
		self.spectral_db = spectral_db
		self.reference_db_conn = sqlite3.connect(reference_db)
		self.spectral_db_conn = sqlite3.connect(spectral_db) if spectral_db else None
		self.read_features()

	def peak_fit(self, t, dt_i, p0="guess"):
		"""
		FeatureAnnotate.peak_fit
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
		FeatureAnnotate.gaussian_fit
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
		return A * np.exp(-(x - B) ** 2 / (2 * C ** 2))

	def fwhm_threshold(self, C, dt_i, fwhm_thresholds=(0.05, 2.5), intensity_threshold=500):
		"""
		FeatureAnnotate.fwhm_thresholds
		description:
			Checks if the full width at half maximum (FWHM) and intensity of the EIM meet
			threshold values.
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

		return fwhm_thresholds[0] < fwhm < fwhm_thresholds[1] and max_intensity > intensity_threshold

	def group_peaks(self, peaks, group_tolerance=1):
		"""
		FeatureAnnotate.group_peaks
		description:
			Groups peaks based on the group tolerance and finds the most intense peak within each group.
			This process is the equivalent of "centroiding" the mass spectra for streamlined analysis.
			It is intended to be used on the reference spectra prior to scoring.
		parameters:
			peaks (list of tuples) -- MS peaks to be grouped.
			group_tolerance (float) -- m/z window used to group peaks. Default is 1.
		returns:
			(list of tuples) -- most intense MS peaks from each group.
		"""
		if not peaks:

			return []

		# Sort peaks by m/z value
		sorted_peaks = sorted(peaks, key=lambda x: x[0])

		# Initialize the first group
		grouped_peaks = [[sorted_peaks[0]]]

		# Iterate over the rest of the sorted peaks to determine if it should be in a new group
		for current_peak in sorted_peaks[1:]:

			# Check if the current peak is within the tolerance of the last group's m/z value
			if abs(current_peak[0] - grouped_peaks[-1][-1][0]) <= group_tolerance:
				grouped_peaks[-1].append(current_peak)
			else:

				# If not, start a new peak group
				grouped_peaks.append([current_peak])

		# For each group, find the peak with the maximum intensity
		most_intense_peaks = [max(group, key=lambda x: x[1]) for group in grouped_peaks]

		return most_intense_peaks

	def gaussian_smooth_pick(self, mz_array, intensity_array, window_len=51, std=7, prominence=0.1):
		"""
		FeatureAnnotate.gaussian_smooth_pick
		description:
			Applies a Gaussian smoothing to the intensity array and identifies peaks.
			This function convolves the intensity array with a Gaussian window to smooth the data.
			It then identifies local maxima (i.e., peaks) in the smoothed data based on a defined prominence.
		parameters:
			mz_array (array-like): array of m/z values corresponding to the intensities.
			intensity_array (array-like): array of intensity values to be smoothed and from which peaks are identified.
			window_len (int) -- length of the gaussian window used for smoothing. Default is 51.
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

	"""def plot_aif(self, mz_array, intensity_array, smoothed_intensity, local_maxima):
		plt.figure(figsize=(10, 6))
		plt.plot(mz_array, intensity_array, label="Raw Extracted AIF", alpha=0.5)
		plt.plot(mz_array, smoothed_intensity, label="Smoothed AIF", alpha=0.7)
		maxima_mz, maxima_intensity = zip(*local_maxima)
		plt.scatter(maxima_mz, maxima_intensity, color="red", marker="x", label="Local Maxima")
		plt.xlabel("m/z")
		plt.ylabel("Intensity")
		plt.title("AIF with Local Maxima")
		plt.legend()
		plt.show()"""

	def weighted_function(self, peak):
		"""
		FeatureAnnotate.weighted_function
		description:
			Calculates the weighted intensity function for cosine similarity scoring of spectral peaks.
			This is the function used by Waters: https://www.nonlinear.com/progenesis/qi/v1.0/faq/database-fragmentation-algorithm.aspx.
			This algorithm also incorporates differential weighting for the fragment ions, similar to MetaboAnnotatoR (Ebbels et al., 2022).
			If the fragment ion in the extracted spectrum is dt-aligned with the precursor ion, it receives the full weight.
			If the fragment ion in the extracted spectrum is from the raw AIF only, it is downweighted by a factor of 2.
		parameters:
			peak (tuple) -- contains the m/z value, intensity, and tag of the peak.
		returns:
			(float) -- Weighted intensity value.
		"""
		mz, intensity, tag = peak
		weight = 1 if tag == "dt" else 0.5 if tag == "mz" else 1

		return (mz * 2) * math.sqrt(intensity) * weight

	def match_mz(self, reference_mz, reference_intensity, extracted_peaks, mz_tolerance=0.025):
		"""
		FeatureAnnotate.match_mz
		description:
			Matches m/z peaks between reference and extracted spectra within a specified tolerance.
			This is the method used by Waters.
		parameters:
			reference_mz (list or numpy.ndarray) -- m/z values from the reference spectrum.
			reference_intensity (list or numpy.ndarray) -- corresponding intensities of the reference spectrum peaks.
			extracted_mz (list) -- m/z value from the extracted spectrum.
			extracted_intensity (list) -- corresponding intensities of the extracted spectrum peaks.
			mz_tolerance (float) -- tolerance for matching peaks. 
		returns:
			(tuple) -- two lists representing the matched intensity vectors from the reference and extracted spectra.
		"""
		# Convert input spectra to lists
		reference_mz_list = reference_mz.tolist() if isinstance(reference_mz, np.ndarray) else list(reference_mz)
		extracted_mz_list = [mz for mz, _, _ in extracted_peaks]

		# Combine and sort all unique m/z values
		all_mz = sorted(set(reference_mz_list + extracted_mz_list))
		"""print("mz list: {}".format(all_mz))"""

		# Initialize lists to hold the matched intensity values
		reference_vector = []
		extracted_vector = []

		# Track matched peaks to prevent rematching
		matched_ref_peaks = set()
		matched_ext_peaks = set()

		for mz in all_mz:

			# Find closest peaks within tolerance, not already matched
			closest_ref_peak = None
			closest_ext_peak = None

			# Check for the closest reference peak within the tolerance that has not been matched yet
			for mz_ref, int_ref in zip(reference_mz, reference_intensity):
				if abs(mz_ref - mz) <= mz_tolerance and mz_ref not in matched_ref_peaks:
					closest_ref_peak = (mz_ref, int_ref)
					break

			# Similar check for the extracted peak
			for mz_ext, int_ext, tag in extracted_peaks:
				if abs(mz_ext - mz) <= mz_tolerance and mz_ext not in matched_ext_peaks:
					closest_ext_peak = (mz_ext, int_ext, tag)
					break

			# Process the matched peaks
			if closest_ref_peak and closest_ext_peak:

				# Add peaks to the matched sets to prevent rematching
				matched_ref_peaks.add(closest_ref_peak[0])
				matched_ext_peaks.add(closest_ext_peak[0])
				
				# Calculate the weighted intensity for each peak
				ref_intensity = self.weighted_function((closest_ref_peak[0], closest_ref_peak[1], "ref"))
				ext_intensity = self.weighted_function(closest_ext_peak)

			# Handle case where there is no matching peak in the extracted spectrum
			# Add a 0 to the extracted vector
			elif closest_ref_peak:
				matched_ref_peaks.add(closest_ref_peak[0])
				ref_intensity = self.weighted_function((closest_ref_peak[0], closest_ref_peak[1], "ref"))
				ext_intensity = 0

			# Handle case where there is no matching peak in the reference spectrum
			# Add a 0 to the reference vector
			elif closest_ext_peak:
				matched_ext_peaks.add(closest_ext_peak[0])
				ref_intensity = 0
				ext_intensity = self.weighted_function(closest_ext_peak)
			else:
				continue

			# Append the calculated intensities to the vectors
			reference_vector.append(ref_intensity)
			extracted_vector.append(ext_intensity)

		return reference_vector, extracted_vector

	def cosine_similarity(self, reference_vector, extracted_vector):
		"""
		FeatureAnnotate.cosine_similarity
		description:
			Calculates the cosine similarity between the reference and extracted spectral vectors.
		parameters:
			reference_vector (list) -- weighted intensity vector of the reference spectrum.
			extracted_vector (list) -- weighted intensity vector of the extracted spectrum.
		returns:
			(float) -- cosine similarity score between the two vectors.
		"""
		norm_ref = np.linalg.norm(reference_vector)
		norm_ext = np.linalg.norm(extracted_vector)
		if norm_ref == 0 or norm_ext == 0:
			return 0
		similarity = (np.dot(reference_vector, extracted_vector) / (norm_ref * norm_ext)) * 100

		return round(similarity, 2)

	def mirror_plot(self, processed_mz, processed_intensity, msms_data, title_mirror, fname_mirror):
		"""
		FeatureAnnotate.mirror_plot
		description:
			Generates a mirror plot displaying the extracted and reference fragmentation spectra.
		parameters:
			processed_mz (list) -- m/z values from the deconvoluted extracted spectrum.
			processed_intensity (list) -- corresponding intensity values from the deconvoluted extracted spectrum.
			msms_data (list of tuples) -- list of m/z and intensity tuples from the reference spectrum.
			title_mirror (str) -- title of the mirror plot.
			fname_mirror (str) -- file name for saving the plot.
		returns:
			(.png) -- image displaying the extracted and reference fragmentation spectra.
		"""
		if len(processed_mz) == 0 or len(msms_data) == 0:
			print("One or both m/z arrays are empty.")

			return

		# Group reference and extracted spectra for ease of viewing
		most_intense_experimental = self.group_peaks(zip(processed_mz, processed_intensity), group_tolerance=1)
		most_intense_reference = self.group_peaks(msms_data, group_tolerance=1)

		# Plot the data
		fig, ax = plt.subplots(figsize=(6.4, 4.8))
		exp_mz, exp_intensity = zip(*most_intense_experimental) if most_intense_experimental else ([], [])
		ref_mz, ref_intensity = zip(*most_intense_reference) if most_intense_reference else ([], [])
		ref_intensity = [-intensity for intensity in ref_intensity]

		stemContainer1 = ax.stem(exp_mz, exp_intensity, linefmt="-b", basefmt=" ", markerfmt=" ", label="Experimental") # Plot experimental MS/MS spectrum
		stemContainer2 = ax.stem(ref_mz, ref_intensity, linefmt="-r", basefmt=" ", markerfmt=" ", label="Reference") # Plot reference MS/MS spectrum
		stemContainer1.stemlines.set_linewidths(1)
		stemContainer2.stemlines.set_linewidths(1) 

		# Annotate the most intense peaks for experimental( i.e, extracted) data
		for mz_val, intensity_val in most_intense_experimental:
			ax.annotate(f"{mz_val:.4f}", xy=(mz_val, intensity_val), xytext=(2, 2),
						textcoords="offset points", ha="center", va="bottom", fontsize=8)

		# Annotate the most intense peaks for reference data
		for mz_val, intensity_val in most_intense_reference:
			ax.annotate(f"{mz_val:.4f}", xy=(mz_val, -intensity_val), xytext=(2, -2),
						textcoords="offset points", ha="center", va="top", fontsize=8)
		legend = ax.legend(loc="best", frameon=True, edgecolor="black", fontsize=10, facecolor="white") # Create legend
		for text in legend.get_texts():
			text.set_fontname("Arial")
		ax.set_title(title_mirror, fontsize=12, fontweight="bold", fontname="Arial")
		ax.set_xlabel("m/z", fontsize=12, fontweight="bold", fontname="Arial")
		ax.set_ylabel("Intensity [%]", fontsize=12, fontweight="bold", fontname="Arial")
		plt.axhline(0, color="black", linewidth=2)
		ax.set_ylim(-110, 110) # Modify y axis range
		max_mz = max(max(exp_mz), max(ref_mz))
		ax.set_xlim(0, max_mz + 50) # Modify x axis range
		current_ticks = plt.gca().get_yticks()
		plt.gca().set_yticks(current_ticks)
		plt.gca().set_yticklabels([abs(tick) for tick in current_ticks])
		ax.spines["top"].set_visible(False)
		ax.spines["right"].set_visible(False)
		ax.spines["left"].set_linewidth(1.5)
		ax.spines["bottom"].set_linewidth(1.5)
		plt.savefig(fname_mirror, dpi=350, bbox_inches="tight")
		plt.close()

	def read_features(self):
		"""
		FeatureAnnotate.read_features
		description:
			Reads Excel (.xlsx) file containing spectral features inputted by the user into a pandas DataFrame.
		"""
		self.feature_df = pd.read_excel(self.feature_list)

	def extract_aif(self, file_name, mobility_function, rt_start, rt_end, dt_start, dt_end):
		"""
		FeatureAnnotate.extract_aif
		description:
			Extracts the raw AIF spectrum for a given spectral feature.
		parameters:
			file_name (str) -- name of the .raw file. 
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
		m, i = rdr.get_spectrum(mobility_function, rt_start, rt_end, dt_min = dt_start, dt_max = dt_end)

		return np.array(m), np.array(i)

	def monoisotopic_peak(self, peaks, target_mz, mz_tolerance=0.025):
		"""
		FeatureAnnotate.monoisotopic_peak
		description:
			Finds the most intense peak in the AIF spectrum that is closest to the theoretical precursor m/z value, within the specified tolerance.
		parameters:
			peaks (list of tuples) -- peaks in the AIF to compare.
			target_mz (float) -- theoretical precursor m/z value. 
			mz_tolerance (float) -- m/z tolerance for considering a peaks as the precursor peak.
		returns:
			(tuple) -- most intense peak that is closest in value to the theoretical precursor.
		"""
		# Filter peaks within the specified m/z tolerance
		potential_precursors = [peak for peak in peaks if abs(peak[0] - target_mz) <= mz_tolerance]

		# If there are no peaks within the tolerance, return None
		if not potential_precursors:
			return None

		# Return the most intense peak within the tolerance range
		return max(potential_precursors, key=lambda x: x[1])

	def deconvolute_spectrum(self, file_name, extracted_mz, extracted_intensity, mz, rt_start, rt_end, dt, mobility_function, ms1_function, reference_spectra, dt_tolerance=0.1, group_tolerance=1):
		"""
		FeatureAnnotate.deconvolute_spectrum
		description:
			Deconvolutes the raw, extracted AIF based on rt and dt alignment.
		parameters:
			file_name (str) -- name of the .raw file.
			extracted_mz (array) -- array of m/z values from the AIF spectrum.
			extracted_intensity (array) -- array of intensity values from the AIF spectrum.
			mz (float) -- precursor m/z value of the spectral feature.
			rt_start (int) -- index corresponding to beginning of EIC peak. 
			rt_end (int) -- index corresponding to end of EIC peak.
			dt (float) -- drift time of the spectral feature.
			mobility_function (int) -- mobility function number.
			ms1_function (int) -- MS1 function number.
			reference_spectra (array-like) -- reference fragmentation m/z and intensity pairs. 
			dt_tolerance (float) -- drift time alignment tolerance for spectral deconvolution of fragment and precursor ions. Default is 0.1 ms. 
			group_tolerance (float) -- tolerance for grouping peaks by m/z in the reference spectrum. Default is 1 Da. 
			
		returns:
			(list) -- deconvoluted AIF spectrum.
		"""
		# Initialize MassLynxReader object
		rdr = MassLynxReader(file_name)
		
		# Group reference spectra for each potential match
		grouped_reference_spectra = {}
		for potential_match in reference_spectra.keys():
			msms_data = reference_spectra[potential_match]
			grouped_reference_peaks = self.group_peaks(msms_data, group_tolerance=1)
			grouped_reference_spectra[potential_match] = grouped_reference_peaks

		processed_peaks = []

		# Smooth and fit the extracted AIF spectrum using Gaussian convolution 
		# This process is analogous to "centroiding" the profile data as outlined in MSnbase
		# The default parmaeters for the smoothing and picking (i.e., window_len=1, std=7, prominence=0.01) pick the most intense peak for each ion distribution.
		# More aggressive smoothing will lead to slightly different m/z peaks being picked.
		identified_peaks, smoothed_intensity = self.gaussian_smooth_pick(extracted_mz, extracted_intensity, window_len=1, std=7, prominence=0.01)

		"""self.plot_aif(extracted_mz, extracted_intensity, smoothed_intensity, identified_peaks)"""

		# Identify the precursor peak in the AIF
		precursor_peak = self.monoisotopic_peak(identified_peaks, mz)
		precursor_dt = dt
		if precursor_peak:

			# Label the precursor peak and ensure it is appended to processed_peaks
			processed_peaks.append((precursor_peak[0], precursor_peak[1], "precursor"))
		print("\n\t\tprecursor peak: {} precursor dt: {}".format(precursor_peak[0], precursor_dt))

		# Iterate over each grouped reference fragment peak to find the closest peak in the AIF
		for _, grouped_ref_peaks in grouped_reference_spectra.items():
			for ref_mz, _ in grouped_ref_peaks:

				# Skip the reference peak closest to the precursor m/z
				if abs(ref_mz - mz) <= 0.025:
					continue

				closest_peak = self.monoisotopic_peak(identified_peaks, ref_mz)	
				if closest_peak:
					mz_val, intensity_val = closest_peak 

					# Skip the identified precursor peak
					if mz_val == precursor_peak[0]:
						continue

					# Check for retention and drift time alignment with the precursor ion
					fragment_dt, fragment_dt_i = rdr.get_filtered_chrom(mobility_function, mz_val, 0.025, rt_start, rt_end)

					# Parameter for smoothing Gaussian function
					t_refined = np.arange(min(fragment_dt), max(fragment_dt), float(0.001))

					# Initialize and fit Gaussian function
					A, B, C = self.peak_fit(fragment_dt, fragment_dt_i)
					fit_i = self.gaussian_fit(t_refined, A, B, C)

					# Apply FWHM and intensity thresholds
					if self.fwhm_threshold(C, fragment_dt_i):
						fwhm = C * 2.255
						"""print("fwhm: {}".format(fwhm))"""
						fitted_fragment_dt = float(round(B, 2))

						print("\n\t\tfragment m/z: {} fragment dt: {}".format(mz_val, fitted_fragment_dt))

						# Apply tags depending on alignment status with the precursor ion
						tag = "dt" if abs(fitted_fragment_dt - precursor_dt) <= dt_tolerance else "mz"
						processed_peaks.append((mz_val, intensity_val, tag))

		# Handle remaining AIF peaks
		# Remove low intensity peaks, default is 5% of precursor peak
		threshold_intensity = precursor_peak[1] * 0.05

		# Remove AIF peaks that have already been processed
		remaining_peaks = [(mz_val, intensity_val) for mz_val, intensity_val in identified_peaks if mz_val <= mz and intensity_val >= threshold_intensity and not any(peak[0] == mz_val for peak in processed_peaks)]
		print(remaining_peaks)

		# Iterate over each extracted peak in the remaining AIF peaks
		grouped_remaining_peaks = {}
		for mz_val, intensity_val in remaining_peaks:
			found_group = False
			for group_mz in grouped_remaining_peaks.keys():
				if abs(mz_val - group_mz) <= 0.025:

					# Group remaining AIF peaks that are within m/z tolerance (i.e., 0.025 Da)
					grouped_remaining_peaks[group_mz].append((mz_val, intensity_val))
					found_group = True 
					break
			if not found_group:
				grouped_remaining_peaks[mz_val] = [(mz_val, intensity_val)]

		print(grouped_remaining_peaks)

		# Process each group of remaining AIF peaks
		for group_mz, group_peaks in grouped_remaining_peaks.items():

			# Choose the most intense peak as the representative peak
			representative_peak = max(group_peaks, key=lambda x: x[1])

			# Attempt to extract precursor retention time-selected EIM
			mz_val, intensity_val = representative_peak
			fragment_dt, fragment_dt_i = rdr.get_filtered_chrom(mobility_function, mz_val, 0.025, rt_start, rt_end)

			# Parameter for smoothing Gaussian curve
			t_refined = np.arange(min(fragment_dt), max(fragment_dt), float(0.001))

			# Initialize and fit Gaussian function
			A, B, C = self.peak_fit(fragment_dt, fragment_dt_i)
			fit_i = self.gaussian_fit(t_refined, A, B, C)

			# Apply FWHM and intensity thresholds
			if self.fwhm_threshold(C, fragment_dt_i) and precursor_dt is not None:
				fitted_fragment_dt = float(round(B, 2))
				print("\n\t\tfragment m/z: {} dt: {}".format(mz_val, fitted_fragment_dt))
				if abs(fitted_fragment_dt - precursor_dt) <= dt_tolerance:

					# If representative peak is drift time-aligned, add all peaks in the group to processed_peaks
					for peak in group_peaks:
						processed_peaks.append((peak[0], peak[1], "dt"))
				
		print("processed peaks: {}".format(processed_peaks))

		return processed_peaks

	def match_features(self, mz_tolerance, rt_tolerance=None, ccs_tolerance=None):
		"""
		FeatureAnnotate.match_features
		description:
			Matches the spectral feature list with potential compounds based on various criteria levels.
		parameters:
			mz_tolerance (float) -- tolerance window around m/z (in Da). 
			rt_tolerance (float) -- tolerance window around rt (optional; read as min if passed).
			ccs_tolerance (float) -- tolerance window around ccs (optional; read as % if passed).
		"""
		# Make sure feature list is populated and read correctly
		if self.feature_df is None:
			raise ValueError("Spectral feature list is empty or not accessible.")

		print("\n...Fetching potential matches from {}...".format(self.reference_db))

		# Create a dictionary to store reference data for each potential match
		match_data = {}

		# Iterate over each row in the pandas DataFrame
		for i, row in self.feature_df.iterrows():
			file_name, mz, gradient, ccs_calibrant, column_type, monoisotopic_mz = row["File Name"], row["Target m/z"], row["Gradient"], row["CCS Calibrant"], row["Column Type"], row["Observed m/z"]
			rt = row.get("Observed Retention Time (min)", None)
			ccs = row.get("Observed CCS (Å²)", None)

			# Print the current spectral feature being processed
			print("\n\traw file: {} m/z: {}".format(file_name, mz))

			# Build the sqlite query based on the criteria
			# In future implementations, this should be made more flexible to accomodate different database formats
			query = "SELECT DISTINCT compound_name, rt, average_ccs FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ?"
			params = [mz - mz_tolerance, mz + mz_tolerance]

			# User has passed values for rt tolerance to FeatureAnnotate.match_features
			if rt_tolerance is not None:

				# Modify the sqlite query to include rt filter
				query += " AND rt BETWEEN ? AND ?"
				params += [rt - rt_tolerance, rt + rt_tolerance]

			# User has passed values for CCS tolerance to FeatureAnnotate.match_features
			if ccs_tolerance is not None:

				# Modify the sqlite query to include CCS filter
				query += " AND average_ccs BETWEEN ? AND ?"
				params += [ccs * (1 - ccs_tolerance), ccs * (1 + ccs_tolerance)]

			# In the current implementation, matches must also have identical experimental conditions
			# In future implementations, this should be made more flexible
			query += "AND gradient = ? AND ccs_calibrant = ? AND column_type = ?"
			params += [gradient, ccs_calibrant, column_type]

			# Initialize sqlite cursor 
			cursor = self.reference_db_conn.cursor()

			# Fetch the results of the query
			cursor.execute(query, params)
			matches = cursor.fetchall()

			# Fetch reference data for each potential match
			detailed_matches = []
			for match in matches:
				compound_name, match_rt, match_ccs = match

				# Calculate the m/z similarity score between potential match m/z (i.e., row["Target m/z"]) and monoisotopic m/z peak of extracted feature
				mz_error = abs(monoisotopic_mz - mz) / mz * 1e6
				mz_similarity_score = min(1 / mz_error, 1) * 100 if mz_error != 0 else 100

				# Calculate error between CCS of extracted spectral feature and reference CCS of potential match
				ccs_difference = abs(match_ccs - ccs) / match_ccs
				ccs_delta = ccs_difference * 100
				ccs_similarity_score = max(1 - ccs_difference, 0) * 100

				# Calculate error between retention time of extracted spectral feature and reference retention time of potential match
				rt_difference = abs(reference_rt - rt) / reference_rt
				rt_delta = rt_difference * 100

				detailed_matches.append({"compound_name": compound_name, "rt": match_rt, "ccs": match_ccs, "mz_error": mz_similarity_score, "ccs_error": ccs_delta, "rt_error": rt_delta, "match_info": f"{compound_name} (Retention Time (min): {match_rt:.2f} | CCS (Å²): {match_ccs:.2f}", "error_info": f"{compound_name} (MASS ERROR (ppm): {mz_error:.2f} | RETENTION TIME ERROR (%): {rt_delta:.2f} | CCS ERROR (%) {ccs_delta:.2f})"})

			# Rank the potential matches based on m/z error
			ranked_matches = sorted(detailed_matches, key=lambda x: abs(x["mz_error"]))

			# Extract the match strings for the DataFrame
			potential_match_strings = [match["match_info"] for match in detailed_matches]
			ranked_match_strings = [match["error_info"] for match in ranked_matches]

			# Extract potential matches from database query	
			self.feature_df.at[i, "Potential Matches"] = ", ".join(potential_match_strings)
			self.feature_df.at[i, "Ranked Matches"] = ", ".join(ranked_match_strings)

		# Filter pandas DataFrame to only include rows with potential matches
		filtered_df = self.feature_df[self.feature_df["Potential Matches"].notna()]

		# Check if filtered pandas DataFrame is empty
		if filtered_df.empty:
			print("No matches found for the given criteria.")
		else:

			# Define the output file name based on matching criteria
			base_name = os.path.splitext(os.path.basename(self.feature_list))[0]
			output_file_suffix = "_mz"
			if rt_tolerance is not None:
				output_file_suffix += "_rt"
			if ccs_tolerance is not None:
				output_file_suffix += "_ccs"
			output_file = "{}{}.xlsx".format(base_name, output_file_suffix)
			filtered_df.to_excel(output_file, index=False)

			print(f"\nSuccessful analysis. View potential matches in {output_file}.\n")

	def normalize_and_process(self, processed_peaks):
		"""
		FeatureAnnotate.normalize_and_process
		description:
			Normalize peaks and remove tags for comparison with reference spectra.
		"""
		if processed_peaks:
			max_intensity = max(intensity for _, intensity, _ in processed_peaks)
			normalized_mz = []
			normalized_intensity = []
			for mz, intensity, tag in processed_peaks:
				normalized_mz.append(mz)
				normalized_intensity.append((intensity / max_intensity) * 100)
			return normalized_mz, normalized_intensity
		else:
			return [], []

	def composite_score(self, mz_similarity_score, cosine_similarity_score, ccs_similarity_score):
		"""
		FeatureAnnotate.composite_score
		description:
			Calculates the composite score for potential matches, reflecting m/z and fragmentation similarity. 
			Method adapted from MetaboAnnotatoR (Ebbels et al., 2022)
		parameters:
			mz_similarity_score (float) -- score for similarity between extracted monoisotopic peak and candidate m/z (i.e., the target m/z).
			cosine_similarity_score (float) -- weighted cosine similarity score that downweights AIF fragment ions relative to dt-aligned fragment ions. 
			ccs_similarity_score (float) -- score for similarity between spectral feature CCS value and reference CCS value for each candidate.
		returns:
			(float) -- composite score. 
		"""
		return mz_similarity_score + cosine_similarity_score + ccs_similarity_score

	def match_msms(self, ms1_function, mobility_function, database_type, mz_tolerance, rt_tolerance, ccs_tolerance):
		"""
		FeatureAnnotate.match_msms
		description:
			Matches MS/MS spectra from extracted spectral features against a reference database
			for features that have been matched on specified m/z, rt, and CCS tolerances. 
			This function is the main method for scoring and ranking potential matches using fragmentation data.
		parameters:
			mobility_function (int) -- mobility function number.
			database_type (str) -- spectral database to use; must be either "experimental" or "theoretical."
			mz_tolerance (float) -- tolerance window around m/z (in Da). 
			rt_tolerance (float) -- tolerance window around rt (read as min if passed).
			ccs_tolerance (float) -- tolerance window around ccs (read as % if passed).
			In future implementations, make the database_type argument more flexible or obsolete. 
		returns:
			(.xlsx) -- Excel spreadsheet containing the potential and ranked matches with corresponding scores.
			(.png) -- image displaying reference and extracted fragmentation spectra for visual inspection. 
		"""
		print("\n...Fetching potential matches and calculating cosine similarity scores using {}...".format(self.spectral_db))

		# Make sure feature list is populated and read correctly
		if self.feature_df is None:
			raise ValueError("Spectral feature list is empty or not accessible.")

		# Validate database type
		if database_type not in ["experimental", "theoretical"]:
			raise ValueError("Invalid database type. Choose 'experimental' or 'theoretical'.")

		# Choose the appropriate database connection based on user input
		# Make this more flexible in future implementations
		db_table = "experimental_msms_data" if database_type == "experimental" else "theoretical_msms_data"

		final_rows = []
		
		# Iterate over each row in the pandas DataFrame
		for i, row in self.feature_df.iterrows():
			file_name, mz, rt, ccs, gradient, ccs_calibrant, column_type, monoisotopic_mz, dt = row["File Name"], row["Target m/z"], row["Observed Retention Time (min)"], row["Observed CCS (Å²)"], row["Gradient"], row["CCS Calibrant"], row["Column Type"], row["Observed m/z"], row["Observed Drift Time (ms)"]
			potential_matches_column = []

			# Print the current peak being analyzed to the terminal
			print("\n\traw file: {} m/z: {} rt: {}".format(file_name, mz, rt))

			# Execute sqlite query based on m/z, rt, and CCS
			cursor = self.reference_db_conn.cursor()
			cursor.execute(f"""SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ? AND rt BETWEEN ? AND ? AND average_ccs BETWEEN ? AND ? AND gradient = ? AND ccs_calibrant = ? AND column_type = ?""", (mz - mz_tolerance, mz + mz_tolerance, rt - rt_tolerance, rt + rt_tolerance, ccs * (1 - ccs_tolerance), ccs * (1 + ccs_tolerance), gradient, ccs_calibrant, column_type))

			# Fetch the results of the query
			matches = cursor.fetchall()

			# Extract potential matches from database query
			if matches:
				potential_matches = [match[0] for match in matches]
				ranked_matches = []
				reference_spectra = {}

				# Extract m/z-selected ion chromatogram (EIC) from the .raw file using dhrmasslynxapi
				rdr = MassLynxReader(file_name)
				rt_array, rt_i_array = rdr.get_chrom(ms1_function, mz, 0.025)
				rt_array = np.array(rt_array)
				rt_i_array = np.array(rt_i_array)

				# Smooth and fit EIC with multigauss module
				peak_indices, _, _, peak_ranges = process_chromatogram(rt_array, rt_i_array, file_name, mz, "sample", generate_images=False)

				# Find the LC peak that is closest to the feature list value
				closest_peak_idx = np.argmin(np.abs(rt_array[peak_indices] - rt))

				# Find the LC peak range (i.e., the start and end indices)
				start_idx, end_idx = peak_ranges[closest_peak_idx]

				# Use the peak indices to define the EIC window
				rt_start = rt_array[start_idx]
				rt_end = rt_array[end_idx]

				# Use the peak indices to extract the (m/z,rt)-selected ion mobilogram (EIM)
				t, dt_i = rdr.get_filtered_chrom(mobility_function, mz, mz_tolerance, rt_min=rt_start, rt_max=rt_end)
				
				# Parameter for smoothing Gaussian function
				t_refined = np.arange(min(t), max(t), float(0.001))

				# Initialize and fit Gaussian function
				A, B, C = self.peak_fit(t, dt_i)
				fit_i = self.gaussian_fit(t_refined, A, B, C)
				fitted_dt = float(round(B, 2))

				# Ensure that the extracted drift time is equal to the drift time of the extracted spectral feature
				# Because this value is extracted using the same functions, it should be identical
				# Using a threshold just in case
				if abs(fitted_dt - dt) <= 0.1:

					# Calculate full width at 1% maximum (FW1M) bounds
					fw1m = 2 * np.sqrt(-2 * C **2 * np.log(0.01))

					# Find the EIM peak range (i.e., the start and end drift times corresponding to the FW1M)
					dt_start = B - (fw1m / 2)
					dt_end = B + (fw1m / 2)

					# Ensure dt_start and dt_end are within the bounds of the original drift time array
					dt_start = max(dt_start, min(t))
					dt_end = min(dt_end, max(t))

				# Iterate over each potential match 
				for potential_match in potential_matches:

					# Fetch reference data for each potential match
					ccs_rt_cursor = self.reference_db_conn.cursor()
					ccs_rt_cursor.execute("SELECT rt, average_ccs FROM qacs_rt_ccs WHERE compound_name = ?", (potential_match,))
					reference_data = ccs_rt_cursor.fetchone()
					if reference_data:
						reference_rt, reference_ccs = reference_data

						# Format the Potential Match output string
						match_info = f"{potential_match} (Retention Time (min): {reference_rt:.2f} | CCS (Å²): {reference_ccs:.2f})"
						potential_matches_column.append(match_info)
					else:
						potential_matches_column.append(potential_match)

					# Fetch reference spectrum for each potential match
					cursor = self.spectral_db_conn.cursor()
					cursor.execute(f"SELECT mz, normalized_intensity FROM {db_table} WHERE compound_name = ?", (potential_match,))

					# Fetch results of the query
					msms_data = cursor.fetchall()

					if msms_data:

						# Store the raw reference fragmentation spectrum for each potential match
						reference_spectra[potential_match] = msms_data
						"""print("reference spectra: {}".format(reference_spectra))"""
					else:
						reference_spectra[potential_match] = []
					if reference_data:
						print("\n\t\t**Potential match found. Performing spectral deconvolution and scoring.**")

						# Calculate error between retention time of extracted spectral feature and reference retention time of potential match
						rt_difference = abs(reference_rt - rt) / reference_rt
						rt_delta = rt_difference * 100

						# Calculate error between CCS of extracted spectral feature and reference CCS of potential match
						ccs_difference = abs(reference_ccs - ccs) / reference_ccs
						ccs_delta = ccs_difference * 100
						ccs_similarity_score = max(1 - ccs_difference, 0) * 100
						"""print("\n\t\tCCS Similarity Score: {}".format(ccs_similarity_score))"""
					else:
						ccs_similarity_score = 0

				# Extrac AIF spectrum for the experimental peak with potential matches
				extracted_mz, extracted_intensity = self.extract_aif(file_name, mobility_function, rt_start, rt_end, dt_start, dt_end)

				# Deconvolute AIF spectrum
				processed_peaks = self.deconvolute_spectrum(file_name, extracted_mz, extracted_intensity, mz, rt_start, rt_end, row.get("Observed Drift Time (ms)", None), mobility_function, ms1_function, reference_spectra, dt_tolerance=0.1, group_tolerance=0.5)

				# Normalize deconvoluted AIF spectrum to most intense peak for fair comparison to normalized reference spectra
				processed_mz, processed_intensity = self.normalize_and_process(processed_peaks)

				# Initialize dictionary to store fragmentation score for each tag
				fragmentation_score = {"precursor": 0, "dt": 0, "mz": 0}
				for _, _, tag in processed_peaks:
					fragmentation_score[tag] += 1

				print("\n\t\tFragmentation Tag Counter: {}".format(fragmentation_score))

				similarity_list = []

				processed_peaks_tuples = [(mz, intensity, tag) for mz, intensity, tag in processed_peaks]

				# Iterate over each potential match 
				for potential_match, msms_data in reference_spectra.items():

					# Construct the vectors for the extracted and reference spectra
					reference_vector, processed_vector = self.match_mz([mz for mz, _ in msms_data], [intensity for _, intensity in msms_data], processed_peaks_tuples, mz_tolerance)
					similarity_score = self.cosine_similarity(reference_vector, processed_vector)

					# Group reference spectra to extract the number of total possible fragments
					grouped_reference_peaks = self.group_peaks(msms_data, group_tolerance=1)

					# Identify the precursor peak in the reference spectrum as the peak closest to the target m/z
					precursor_peak = self.monoisotopic_peak(grouped_reference_peaks, mz)

					# Count the total possible fragments in the reference spectrum (excluding the precursor peak)
					total_possible_fragments = len([peak for peak in grouped_reference_peaks if peak[0] != precursor_peak[0]])

					# Count the number of matched fragment peaks (excluding the precursor peak)
					matched_fragment_count = sum(1 for _, _, tag in processed_peaks_tuples if tag in ["mz", "dt"])

					# Calculate m/z similarity score between potential match m/z (i.e, row["Target m/z"]) and monoisotopic m/z peak of extracted feature
					mz_error = abs(monoisotopic_mz - mz) / mz * 1e6 
					mz_similarity_score = min(1 / mz_error, 1) * 100 if mz_error != 0 else 100

					# Calculate composite score
					composite_score = self.composite_score(mz_similarity_score, similarity_score, ccs_similarity_score)

					# Store the composite score to show both the raw score and its percentage out of the theoretical maximum score of 300
					formatted_composite_score = f"{composite_score:.2f} / {(composite_score / 300) * 100:.2f}%"
					formatted_ranked_match = f"{potential_match} (COMPOSITE: {formatted_composite_score} | COSINE: {similarity_score:.2f} (MATCHED FRAGMENTS: {matched_fragment_count} / {total_possible_fragments}) | MASS ERROR (ppm): {mz_error:.2f} | RETENTION TIME ERROR (%): {rt_delta:.2f} | CCS ERROR (%): {ccs_delta:.2f})"

					similarity_list.append((composite_score, similarity_score, formatted_ranked_match, potential_match))

					# Create "Fragmentation Spectra" folder if not already present
					spectra_directory = "Fragmentation Spectra"
					if not os.path.exists(spectra_directory):
						os.makedirs(spectra_directory)

					# Generate mirror plot
					title_mirror = "Experimental vs. Reference MS/MS Spectra \nPotential Match: {} (Score: {}) ".format(potential_match, "0" if np.isnan(similarity_score) else "{:.2f}".format(similarity_score))
					if database_type == "experimental":
						fname_mirror = "{}/{}_{}_{}_{:.2f}_Experimental_MSMS.png".format(spectra_directory, file_name, mz, potential_match, rt)
						self.mirror_plot(list(processed_mz), list(processed_intensity), msms_data, title_mirror, fname_mirror)
					elif database_type == "theoretical":
						fname_mirror = "{}/{}_{}_{}_{:.2f}_Theoretical_MSMS.png".format(spectra_directory, file_name, mz, potential_match, rt)
						self.mirror_plot(list(processed_mz), list(processed_intensity), msms_data, title_mirror, fname_mirror)

				# Sort and format the potential matches by composite score
				formatted_ranked_matches = ", ".join([t[2] for t in sorted(similarity_list, key=lambda x: x[0], reverse=True)]) if similarity_list else ""
				formatted_potential_matches = ", ".join(potential_matches_column) if potential_matches_column else None

				# Update the pandas DataFrame
				self.feature_df.at[i, "Potential Matches"] = formatted_potential_matches
				self.feature_df.at[i, "Ranked Matches"] = formatted_ranked_matches
				"""print(self.feature_df[["potential_matches", "ranked_matches"]].head())"""

		# Filter pandas DataFrame to remove rows with no potential matches or ranked matches
		self.feature_df.replace("nan", np.nan, inplace=True)
		filtered_df = self.feature_df.dropna(subset=["Potential Matches", "Ranked Matches"])
		filtered_df = filtered_df[filtered_df["Potential Matches"] != ""]
		filtered_df = filtered_df[filtered_df["Ranked Matches"] != ""]
		"""print(filtered_df[["potential_matches", "ranked_matches"]].head())"""

		# Choose output file name based on database type
		base_name = os.path.splitext(os.path.basename(self.feature_list))[0]
		output_file_suffix = "_matches_ranked_experimental.xlsx" if database_type == "experimental" else "_matched_ranked_theoretical.xlsx"
		output_file = base_name + output_file_suffix	

		# Export the pandas DataFrame to an Excel (.xlsx) file
		filtered_df.to_excel(output_file, index=False)

		print(f"\nSuccessful analysis. View ranked matches in {output_file}.\n")

	def __del__(self):
		# Close the database connections
		self.reference_db_conn.close()
		if self.spectral_db_conn:
			self.spectral_db_conn.close()