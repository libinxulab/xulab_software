import pandas as pd 
import sqlite3
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.optimize import curve_fit
from multigauss import process_data
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
			This process is the equivalent of de-isotoping the mass spectra for streamlined analysis.
		parameters:
			peaks (list of tuples) -- MS peaks to be grouped.
			group_tolerance (float) -- m/z window used to group peaks.
		returns:
			(list of tuple) -- most intense MS peaks from each group.
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
			reference_mz (list) -- m/z values from the reference spectrum.
			reference_intensity (list) -- corresponding intensities of the reference spectrum peaks.
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

	def mirror_plot(self, processed_mz, processed_intensity, reference_peaks, title_mirror, fname_mirror):
		"""
		FeatureAnnotate.mirror_plot
		description:
			Generates a mirror plot displaying the extracted and reference spectra.
		parameters:
			processed_mz (list) -- m/z values from the deconvoluted extracted spectrum.
			processed_intensity (list) -- corresponding intensity values from the deconvoluted extracted spectrum.
			reference_peaks (list of tuples) -- list of m/z and intensity tuples from the reference spectrum.
			title_mirror (str) -- title of the mirror plot.
			fname_mirror (str) -- file name for saving the plot.
		"""
		if len(processed_mz) == 0 or len(reference_peaks) == 0:
			print("One or both m/z arrays are empty.")
			return
		most_intense_experimental = self.group_peaks(zip(processed_mz, processed_intensity), group_tolerance=1)
		most_intense_reference = self.group_peaks(reference_peaks, group_tolerance=1)

		# Plotting
		fig, ax = plt.subplots(figsize=(6.4, 4.8))
		exp_mz, exp_intensity = zip(*most_intense_experimental) if most_intense_experimental else ([], [])
		ref_mz, ref_intensity = zip(*most_intense_reference) if most_intense_reference else ([], [])
		ref_intensity = [-intensity for intensity in ref_intensity]

		stemContainer1 = ax.stem(exp_mz, exp_intensity, linefmt="-b", basefmt=" ", markerfmt=" ", label="Experimental") # Plot experimental MS/MS spectrum
		stemContainer2 = ax.stem(ref_mz, ref_intensity, linefmt="-r", basefmt=" ", markerfmt=" ", label="Reference") # Plot reference MS/MS spectrum
		stemContainer1.stemlines.set_linewidths(1)
		stemContainer2.stemlines.set_linewidths(1) 

		# Annotate the most intense peaks for experimental data
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

	def extract_aif(self, file_name, mobility_function, mz, rt, dt):
		"""
		FeatureAnnotate.extract_aif
		description:
			Extracts the raw AIF spectrum for a given spectral feature.
		parameters:
			file_name (str) -- name of the raw data file.
			mz (float) -- m/z value of the spectral feature used to extract the AIF.
			rt (float) -- retention time of the extracted spectral feature.
			dt (float) -- drift time of the extracted spectral feature, if available.
			mobility_function (int) -- mobility function number.
		returns:
			(tuple) -- arrays of m/z values and corresponding intensities.
		"""
		# Initialize MassLynxReader object
		rdr = MassLynxReader(file_name)
		if not pd.isna(dt):
			dt = float(dt)

			# Extract raw (rt,dt)-selected AIF from the spectral feature
			m, i = rdr.get_spectrum(mobility_function, rt - 0.01, rt + 0.01, dt_min = dt - 0.1, dt_max = dt + 0.1)
		
		# Extract rt-selected AIF from the spectral feature if there is no dt available
		else:
			m, i = rdr.get_spectrum(mobility_function - 1, rt - 0.01, rt + 0.01)

		return np.array(m), np.array(i)

	def find_closest_peak(self, peaks, target_mz):
		"""
		FeatureAnnotate.find_closest_peak
		description:
			Finds the peak in the AIF spectrum that is closest to the theoretical precursor m/z value.
		parameters:
			peaks (list of tuples) -- peaks in the AIF to compare.
			target_mz (float) -- theoretical precursor m/z value. 
		returns:
			(tuple) -- peak that is closest in value to the theoretical precursor.
		"""
		# Return None if the peak list is empty
		if not peaks:
			return None

		return min(peaks, key=lambda x: abs(x[0] - target_mz))

	def deconvolute_spectrum(self, file_name, extracted_mz, extracted_intensity, mz, rt, dt, mobility_function, ms1_function, reference_spectra, dt_tolerance=0.2, group_tolerance=1):
		"""
		FeatureAnnotate.deconvolute_spectrum
		description:
			Deconvolutes the raw, extracted AIF based on rt and dt alignment.
		parameters:
			extracted_mz (array) -- array of m/z values from the AIF spectrum.
			extracted_intensity (array) -- array of intensity values from the AIF spectrum.
			mz (float) -- precursor m/z value of the spectral feature.
			dt (float) -- drift time of the spectral feature.
			dt_tolerance (float) -- drift time alignment tolerance for spectral deconvolution.
			group_tolerance (float) -- tolerance for grouping peaks by m/z in the AIF.
			mobility_function (int) -- mobility function number.
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

		# Identify the precursor peak in the AIF
		precursor_peak = self.find_closest_peak(zip(extracted_mz, extracted_intensity), mz)
		precursor_dt = dt
		if precursor_peak:
			processed_peaks.append((precursor_peak[0], precursor_peak[1], "precursor"))
		print("\n\t\tprecursor peak: {} precursor dt: {}".format(precursor_peak[0], precursor_dt))

		# Iterate over each grouped reference fragment peak to find the closest peak in the AIF
		for _, grouped_ref_peaks in grouped_reference_spectra.items():
			for ref_mz, _ in grouped_ref_peaks:

				# Skip the reference peak closest to the precursor m/z
				if abs(ref_mz - mz) <= 0.025:
					continue
				closest_peak = self.find_closest_peak(zip(extracted_mz, extracted_intensity), ref_mz)	
				if closest_peak:
					mz_val, intensity_val = closest_peak 

					# Skip the identified precursor peak
					if mz_val == precursor_peak[0]:
						continue

					# Check for retention and drift time alignment with the precursor ion
					# In future implementations, link retention time bounds to identified EIC indices
					fragment_dt, fragment_dt_i = rdr.get_filtered_chrom(mobility_function, mz_val, 0.025, rt_min=rt-0.05, rt_max=rt+0.05)

					# Parameter for smoothing Gaussian function
					t_refined = np.arange(min(fragment_dt), max(fragment_dt), float(0.001))

					# Initialize and fit Gaussian function
					A, B, C = self.peak_fit(fragment_dt, fragment_dt_i)
					fit_i = self.gaussian_fit(t_refined, A, B, C)

					# Apply FWHM and intensity thresholds
					if self.fwhm_threshold(C, fragment_dt_i):
						fitted_fragment_dt = float(round(B, 2))

						print("\n\t\tfragment m/z: {} fragment dt: {}".format(mz_val, fitted_fragment_dt))

						# Apply tags depending on alignment status
						tag = "dt" if abs(fitted_fragment_dt - precursor_dt) <= dt_tolerance else "mz"
						processed_peaks.append((mz_val, intensity_val, tag))

		# Group remaining AIF fragment peaks
		remaining_peaks = [(mz_val, intensity_val) for mz_val, intensity_val in zip(extracted_mz, extracted_intensity) if mz_val <= mz and not any(peak[0] == mz_val for peak in processed_peaks)]
		grouped_remaining_peaks = self.group_peaks(remaining_peaks, group_tolerance)
		"""print("grouped remaining peaks: {}".format(grouped_remaining_peaks))"""

		# Iterate over each extracted peak in the remaining grouped AIF spectrum
		for mz_val, intensity_val in grouped_remaining_peaks:

			# Attempt to extract precursor retention time-selected EIM
			# In future implementations, link retention time bounds to EIC indices
			fragment_dt, fragment_dt_i = rdr.get_filtered_chrom(mobility_function, mz_val, 0.025, rt_min=rt-0.05, rt_max=rt+0.05)

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
					processed_peaks.append((mz_val, intensity_val, "dt"))
			
		"""print("processed peaks: {}".format(processed_peaks))"""
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

		# Iterate over each row in the pandas DataFrame
		for i, row in self.feature_df.iterrows():
			file_name, mz, gradient, ccs_calibrant, column_type = row["file_name"], row["mz"], row["gradient"], row["ccs_calibrant"], row["column_type"]
			rt = row.get("rt", None)
			ccs = row.get("ccs", None)

			# Print the current spectral feature being processed
			print("\n\traw file: {} m/z: {}".format(file_name, mz))

			# Build the sqlite quuery based on desired criteria
			query = "SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ?"
			params = [mz - mz_tolerance, mz + mz_tolerance]

			# Build the sqlite query based on the criteria
			# In future implementations, this should be made more flexible to accomodate different database formats
			query = "SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ?"
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

			# Extract potential matches from database query
			if matches:
				self.feature_df.at[i, "potential_matches"] = ", ".join(set(match[0] for match in matches))
			else:
				self.feature_df.at[i, "potential_matches"] = None

		# Filter pandas DataFrame to only include rows with potential matches
		filtered_df = self.feature_df[self.feature_df["potential_matches"].notna()]

		# Check if filtered pandas DataFrame is empty
		if filtered_df.empty:
			print("No matches found for the given criteria.")
		else:

			# Define the output file name based on matching criteria
			base_name = os.path.splitext(os.path.basename(self.feature_list))[0]
			if rt_tolerance is None and ccs_tolerance is None:
				output_file = base_name + "_matches_mz.xlsx"
			elif rt_tolerance is not None and ccs_tolerance is None:
				output_file = base_name + "_matches_mz_rt.xlsx"
			elif rt_tolerance is None and ccs_tolerance is not None:
				output_file = base_name + "_matches_mz_ccs.xlsx"
			elif rt_tolerance is not None and ccs_tolerance is not None:
				output_file = self.feature_list[:-14] + "_matches_mz_rt_ccs.xlsx"
			filtered_df.to_excel(output_file, index=False)

			print(f"\nSuccessful analysis. View potential matches in {output_file}.\n")

	def reference_match(self, mz_val, reference_spectra, mz_tolerance=0.025):
		"""
		FeatureAnnotate.reference_match
		description:
			Checks if a given extracted fragment m/z value has a corresponding peak in the reference spectrum within a specified m/z tolerance.
		parameters:
			mz_val (float) -- the m/z value to be matched.
			reference_spectra (dict) -- dictionary where keys are potential compound names and values are lists of tuples representing m/z, intensity peaks 
			mz_tolerance (float) -- the m/z tolerance for matching peaks. Default is 0.025.
		returns:
			(bool) -- True if a matching peak is found in the reference spectrum, False otherwise. 
		"""
		for _, msms_data in reference_spectra.items():
			for ref_mz, _ in msms_data:
				if abs(ref_mz - mz_val) <= mz_tolerance:
					return True 
		return False

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

	def composite_score(self, mz_similarity_score, fragmentation_score, ccs_similarity_score):
		"""
		FeatureAnnotate.composite_score
		description:
			Calculates the composite score for potential matches, reflecting m/z and fragmentation similarity. 
			Method adapted from MetaboAnnotatoR (Ebbels et al., 2022)
		parameters:
			mz_similarity_score (float) -- score for similarity between extracted monoisotopic peak and candidate m/z (i.e., the target m/z).
			fragmentation_score (float) -- weighted cosine similarity score that downweights AIF fragment ions relative to dt-aligned fragment ions. 
			dt_tag_weight (float) -- weighting factor for dt-aligned fragment ions.
			mz_tag_weight (float) -- weighting factor for non-dt aligned fragment ions that were pulled from AIF. 
		returns:
			(float) -- composite score. 
		"""
		dt_tag_weight = 1
		mz_tag_weight = 0.5

		dt_aligned_score = sum(fragmentation_score.get(tag, 0) * dt_tag_weight for tag in fragmentation_score if tag == "dt")
		mz_aligned_score = sum(fragmentation_score.get(tag, 0) * mz_tag_weight for tag in fragmentation_score if tag == "mz")

		# Combined the scores
		total_fragmentation_score = dt_aligned_score + mz_aligned_score

		return mz_similarity_score + total_fragmentation_score + ccs_similarity_score

	def match_msms(self, ms1_function, mobility_function, database_type, mz_tolerance, rt_tolerance, ccs_tolerance):
		"""
		FeatureAnnotate.match_msms
		description:
			Matches MS/MS spectra from extracted spectral features against a reference database
			for features that have been matched on specified m/z, rt, and CCS tolerances.
		parameters:
			mobility_function (int) -- mobility function number.
			database_type (str) -- spectral database to use; must be either "experimental" or "theoretical."
			mz_tolerance (float) -- tolerance window around m/z (in Da). 
			rt_tolerance (float) -- tolerance window around rt (read as min if passed).
			ccs_tolerance (float) -- tolerance window around ccs (read as % if passed).
			In future implementations, make the database_type argument more flexible or obsolete. 
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
			file_name, mz, rt, ccs, gradient, ccs_calibrant, column_type, monoisotopic_mz = row["file_name"], row["mz"], row["rt"], row["ccs"], row["gradient"], row["ccs_calibrant"], row["column_type"], row["monoisotopic_mz"]
			potential_matches_column = []

			# Print the current peak being analyzed to the terminal
			print("\n\traw file:{} m/z: {} rt: {}".format(file_name, mz, rt))

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

				# Fetch reference spectrum for each potential match
				for potential_match in potential_matches:
					cursor = self.spectral_db_conn.cursor()
					cursor.execute(f"SELECT mz, normalized_intensity FROM {db_table} WHERE compound_name = ?", (potential_match,))

					# Fetch results of the query
					msms_data = cursor.fetchall()

					if msms_data:

						# Group and find the most intense reference peaks for fair comparison to extracted peaks
						grouped_reference_peaks = self.group_peaks(msms_data, group_tolerance=1)
						reference_spectra[potential_match] = grouped_reference_peaks
						"""print("reference spectra: {}".format(reference_spectra))
						print("grouped reference spectra: {}".format(grouped_reference_peaks))"""
					else:
						reference_mz, reference_intensity = [], []

					# Fetch reference CCS value for each potential match
					ccs_cursor = self.reference_db_conn.cursor()
					ccs_cursor.execute("SELECT average_ccs FROM qacs_rt_ccs WHERE compound_name = ?", (potential_match,))
					reference_ccs_data = ccs_cursor.fetchone()

					if reference_ccs_data:
						reference_ccs = reference_ccs_data[0]
						"""print("reference_ccs: {}".format(reference_ccs))"""

						# Calculate error between CCS of extracted spectral feature and reference CCS of potential match
						ccs_difference = abs(reference_ccs - ccs) / reference_ccs
						ccs_similarity_score = max(1 - ccs_difference, 0)
						print("\n\t\t**Potential match found. Performing spectral deconvolution and scoring.**")
						print("\n\t\tCCS Similarity Score: {}".format(ccs_similarity_score))
					else:
						ccs_similarity_score = 0

				# Extrac AIF spectrum for the experimental peak with potential matches
				extracted_mz, extracted_intensity = self.extract_aif(file_name, mobility_function, mz, rt, row.get("dt", None))

				# Deconvolute AIF spectrum
				processed_peaks = self.deconvolute_spectrum(file_name, extracted_mz, extracted_intensity, mz, rt, row.get("dt", None), mobility_function, ms1_function, reference_spectra, dt_tolerance=0.1, group_tolerance=1)

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

					# Calculate m/z similarity score between potential match m/z (i.e, row["mz"]) and monoisotopic m/z peak of extracted feature
					mz_error = abs(monoisotopic_mz - mz) / mz * 1e6 
					mz_similarity_score = min(1 / mz_error, 1) if mz_error != 0 else 1

					"""print("mz similarity score: {}".format(mz_similarity_score))"""

					# Calculate composite score
					composite_score = self.composite_score(mz_similarity_score, fragmentation_score, ccs_similarity_score)

					similarity_list.append((composite_score, similarity_score, potential_match))

					# Create "Fragmentation Spectra" folder if not already present
					spectra_directory = "Fragmentation Spectra"
					if not os.path.exists(spectra_directory):
						os.makedirs(spectra_directory)

					# Generate mirror plot
					title_mirror = "Experimental vs. Reference MS/MS Spectra \nPotential Match: {} (Score: {}) ".format(potential_match, "0" if np.isnan(similarity_score) else "{:.2f}".format(similarity_score))
					if database_type == "experimental":
						fname_mirror = "{}/{}_{}_{}_{:.2f}_Experimental_MSMS.png".format(spectra_directory, file_name, mz, potential_match, rt)
						self.mirror_plot(list(processed_mz), list(processed_intensity), grouped_reference_peaks, title_mirror, fname_mirror)
					elif database_type == "theoretical":
						fname_mirror = "{}/{}_{}_{}_{:.2f}_Theoretical_MSMS.png".format(spectra_directory, file_name, mz, potential_match, rt)
						self.mirror_plot(list(processed_mz), list(processed_intensity), grouped_reference_peaks, title_mirror, fname_mirror)

				# Sort and format the potential matches
				formatted_ranked_matches = ", ".join([f"{t[2]} (CS: {t[0]:.2f}, SIM: {t[1]:.2f})" for t in sorted(similarity_list, key=lambda x: x[0], reverse=True)]) if similarity_list else ""
				formatted_potential_matches = ", ".join(potential_matches) if potential_matches else None

				# Update the pandas DataFrame
				self.feature_df.at[i, "potential_matches"] = formatted_potential_matches
				self.feature_df.at[i, "ranked_matches"] = formatted_ranked_matches
				"""print(self.feature_df[["potential_matches", "ranked_matches"]].head())"""

		# Filter pandas DataFrame to remove rows with no potential matches or ranked matches
		self.feature_df.replace("nan", np.nan, inplace=True)
		filtered_df = self.feature_df.dropna(subset=["potential_matches", "ranked_matches"])
		filtered_df = filtered_df[filtered_df["potential_matches"] != ""]
		filtered_df = filtered_df[filtered_df["ranked_matches"] != ""]
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