import pandas as pd
import sqlite3
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from dhrmasslynxapi.reader import MassLynxReader

# Create "Extracted Mobilograms" folder if it does not exist
spectra_directory = "Fragmentation Spectra"
if not os.path.exists(spectra_directory):
    os.makedirs(spectra_directory)

# Enter path to the Excel spreadsheet containing processed and filtered data
input_file = "2023_07_27_RO1_pooled_feces_30B_90B_filtered.xlsx"

# Mapping numbers to matching criteria for user input
matching_criteria_dict = {
    "1": "mz",
    "2": "mz_rt",
    "3": "mz_ccs",
    "4": "mz_rt_ccs",
    "5": "mz_rt_ccs_msms"
}

# Prompt user for desired matching criteria
print("\nPlease enter the number corresponding to the desired matching criteria:")
print("\n\t> 1: m/z")
print("\n\t> 2: m/z and rt")
print("\n\t> 3: m/z and ccs")
print("\n\t> 4: m/z, rt, and ccs")
print("\n\t> 5: m/z, rt, ccs, and MS/MS similarity score")
selection = input("> ")

# Validate the user input
while selection not in matching_criteria_dict:
    print("Invalid selection. Please enter a number between 1 and 5.")
    print("\n")
    selection = input(">\t")

matching_criteria = matching_criteria_dict[selection]

# If user selects "5: m/z, rt, ccs, and MS/MS similarity score," ask which reference database to use
if matching_criteria == "mz_rt_ccs_msms":
    database_type = input("\nPlease select the desired reference MS/MS database:\n\n\t> 1: Experimental\n\n\t> 2: Theoretical\n\n> ")
    if database_type == "1":
        database_type = "qacs_experimental_msms"
    elif database_type == "2":
        database_type = "qacs_theoretical_msms"

# Read input file into a pandas dataframe
df_input = pd.read_excel(input_file)

# Connect to the qacs database
conn = sqlite3.connect("qacs.db")
c = conn.cursor()

# Connect to the qacs_experimental_msms database
conn_msms = sqlite3.connect("qacs_experimental_msms.db")
c_msms = conn_msms.cursor()

# Connect to the qacs_theoretical_msms database
conn_theoretical_msms = sqlite3.connect("qacs_theoretical_msms.db")
c_theoretical_msms = conn_theoretical_msms.cursor()

# Filter peaks based on their intensity relative to the base peak 
def filter_peaks_by_intensity(peaks, base_peak_intensity, threshold_percentage):
    threshold_intensity = base_peak_intensity * threshold_percentage / 100
    return [(mz, intensity) for mz, intensity in peaks if intensity >= threshold_intensity]

# Calculate the closeness of the m/z match within the tolerance window
def calculate_closeness(extracted_mz, reference_mz, tolerance):
    difference = abs(extracted_mz - reference_mz)
    # Calculate closeness as a ratio of the difference to the tolerance
    closeness = max(0, 1 - difference / tolerance) 
    return closeness 

# Calculate the weighted intensity using the closeness of the m/z match
def calculate_weighted_intensity_function(mz, intensity, closeness):
    weighting_factor = 0.5 + 0.5 * closeness * 2
    return (mz) * (intensity) * weighting_factor # Adjust contributions from mz and intensity 

# Set m/z tolerance for peak alignment and matching
mz_tolerance = 0.025

# Set group tolerance for dealing with closely spaced peaks
group_tolerance = 25

# Match extracted m/z peaks to the reference spectrum
# Peaks that fall outside of the tolerance range do not contribute to the score
def match_mz(reference_peaks, extracted_peaks, mz_tolerance, group_tolerance):
    # First, filter out the extracted peaks that are within the mz_tolerance of any reference peak
    tolerance_filtered_extracted_peaks = [
        (mz, intensity) for mz, intensity in extracted_peaks
        if any(abs(mz - ref_mz) <= mz_tolerance for ref_mz, _ in reference_peaks)
    ]
    # Group peaks and select the most intense peak within each group for both reference and extracted peaks
    grouped_reference_peaks = group_and_find_most_intense_peaks(reference_peaks, group_tolerance)
    grouped_extracted_peaks = group_and_find_most_intense_peaks(tolerance_filtered_extracted_peaks, group_tolerance)
    reference_vector = []
    extracted_vector = []
    # Iterate over the grouped reference peaks
    for ref_mz, ref_intensity in grouped_reference_peaks:
        # Find the closest extracted peak within the mz tolerance
        closest_peak = min(
            (peak for peak in grouped_extracted_peaks if abs(ref_mz - peak[0]) <= mz_tolerance),
            key=lambda x: abs(ref_mz - x[0]),
            default=(None, None)
        )
        closeness = calculate_closeness(ref_mz, closest_peak[0], mz_tolerance) if closest_peak[0] is not None else 0
        if closest_peak[0] is not None:
            # Calculate weighted intensities
            ref_weighted_intensity = calculate_weighted_intensity_function(ref_mz, ref_intensity, closeness)
            ext_weighted_intensity = calculate_weighted_intensity_function(closest_peak[0], closest_peak[1], closeness)
            extracted_vector.append(ext_weighted_intensity)
        else:
            # No peak close enough in the extracted spectrum, so use zero intensity
            extracted_vector.append(0)
            ref_weighted_intensity = calculate_weighted_intensity_function(ref_mz, ref_intensity, closeness)      
        reference_vector.append(ref_weighted_intensity)        
    # Calculate the similarity score
    return similarity_score(reference_vector, extracted_vector)

# Remove ions that are greater than 2 Da above the molecular ion (i.e. "mz")
def filter_ions_by_mz(peaks, molecular_ion_mz):
    return [(mz, intensity) for mz, intensity in peaks if mz <= molecular_ion_mz + 2]

# Calculate the cosine similarity score between the two spectra
def similarity_score(reference_vector, extracted_vector):
    norm_ref = np.linalg.norm(reference_vector)
    norm_ext = np.linalg.norm(extracted_vector)  
    # Normalize vectors to unit vectors
    normalized_ref_vector = reference_vector / norm_ref if norm_ref != 0 else reference_vector
    normalized_ext_vector = extracted_vector / norm_ext if norm_ext != 0 else extracted_vector
    # Calculate cosine similarity score
    similarity = np.dot(normalized_ref_vector, normalized_ext_vector)
    return similarity * 100

# Set global font conditions for figures
params = {"font.family": "Arial",
          "font.weight": "bold"}
plt.rcParams.update(params)

def group_and_find_most_intense_peaks(peaks, group_tolerance):
    if not peaks:
        return []
    # Sort peaks by m/z value
    sorted_peaks = sorted(peaks, key=lambda x: x[0])
    # Initialize the first group
    grouped_peaks = [[sorted_peaks[0]]]
    for current_peak in sorted_peaks[1:]:
        # Check if the current peak is within the tolerance of the last group's m/z value
        if abs(current_peak[0] - grouped_peaks[-1][-1][0]) <= group_tolerance:
            grouped_peaks[-1].append(current_peak)
        else:
            # If not, start a new group
            grouped_peaks.append([current_peak])
    # For each group, find the peak with the maximum intensity
    most_intense_peaks = [max(group, key=lambda x: x[1]) for group in grouped_peaks]
    return most_intense_peaks

def mirror_plot(filtered_mz, filtered_intensity, reference_peaks, title_mirror, fname_mirror, mz_tolerance, group_tolerance):
    if len(filtered_mz) == 0 or len(reference_peaks) == 0:
        print("One or both m/z arrays are empty.")
        return
    reference_mz, reference_intensity = zip(*reference_peaks)
    reference_intensity = [-intensity for intensity in reference_intensity]
    tolerance_filtered_peaks = [(mz, intensity) for mz, intensity in zip(filtered_mz, filtered_intensity)
                                if any(abs(mz - ref_mz) <= mz_tolerance for ref_mz in reference_mz)]
    if not tolerance_filtered_peaks:
        print("No extracted peaks within the m/z tolerance range.")
        return
    tolerance_filtered_mz, tolerance_filtered_intensity = zip(*tolerance_filtered_peaks)
    # Group peaks and select the most intense peak within each group
    most_intense_experimental = group_and_find_most_intense_peaks(tolerance_filtered_peaks, group_tolerance)
    most_intense_reference = group_and_find_most_intense_peaks(reference_peaks, group_tolerance)
    # Plotting
    fig, ax = plt.subplots(figsize=(6.4, 4.8))
    stemContainer1 = ax.stem(tolerance_filtered_mz, tolerance_filtered_intensity, linefmt="-b", basefmt=" ", markerfmt=" ", label="Experimental") # Plot experimental MS/MS spectrum
    stemContainer2 = ax.stem(reference_mz, reference_intensity, linefmt="-r", basefmt=" ", markerfmt=" ", label="Reference") # Plot reference MS/MS spectrum
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
    max_mz = max(max(filtered_mz), max(reference_mz))
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

final_rows = []

if matching_criteria in ["mz", "mz_rt", "mz_ccs", "mz_rt_ccs"]: 

    print("\n...Fetching potential matches...")

    rows_with_matches = []

    # Iterate over each row in the Excel file
    for i, row in df_input.iterrows():

        # Extract the values for m/z, rt, ccs, ccs_calibrant, gradient, and column_type for SQLite3 query
        # Extract the remaining metadata to populate dataframe
        file_name = row["file_name"]
        sample_type = row["sample_type"]
        mz = row["mz"]
        ccs_calibrant = row["ccs_calibrant"]
        gradient = row["gradient"]
        column_type = row["column_type"]
        dt = row["dt"]
        ccs = row["ccs"]
        rt = row["rt"]
        peak_area = row["peak_area"]

        # Print the current peak being analyzed to the terminal
        print(f"\n\traw file: {file_name} m/z: {mz}")

        # Query the qacs database for matches based on selected matching criteria
        # Default criteria are: +/- 0.025 Da for m/z,  +/- 0.2s for rt, +/- 3% for ccs
        # The match must also have identical experimental conditions (i.e. gradient, column type, CCS calibrant)
        if matching_criteria == "mz":
            c.execute('''SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ?
                         AND gradient = ? AND ccs_calibrant = ? AND column_type = ?''', 
                     (mz - 0.025, mz + 0.025, gradient, ccs_calibrant, column_type))
        elif matching_criteria == "mz_rt":
            c.execute('''SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ?
                         AND rt BETWEEN ? AND ? AND gradient = ? AND ccs_calibrant = ?
                         AND column_type = ?''', (mz - 0.025, mz + 0.025, rt - 0.2, rt + 0.2, 
                         gradient, ccs_calibrant, column_type))
        elif matching_criteria == "mz_ccs":
            c.execute('''SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ?
                         AND average_ccs BETWEEN ? AND ? AND gradient = ? AND 
                         ccs_calibrant = ? AND column_type = ?''', (mz - 0.025, mz + 0.025, ccs * 0.97,
                         ccs * 1.03, gradient, ccs_calibrant, column_type))
        elif matching_criteria == "mz_rt_ccs":
            c.execute('''SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ?
                         AND rt BETWEEN ? AND ? AND average_ccs BETWEEN ? AND ?
                         AND gradient = ? AND ccs_calibrant = ? AND column_type = ?''', 
                         (mz - 0.025, mz + 0.025, rt - 0.2, rt + 0.2, ccs * 0.97, ccs * 1.03,
                          gradient, ccs_calibrant, column_type))

        # Get the result of the query
        matches = c.fetchall()

        # Extract potential matches from database query
        if matches:
            row["potential_matches"] = " ,  ".join(set([match[0] for match in matches]))
            rows_with_matches.append(row)

    # Save updated input dataframe to a new Excel file
    if matching_criteria == "mz":
        output_file = input_file[:-14] + "_matches_mz.xlsx"
    elif matching_criteria == "mz_rt":
        output_file = input_file[:-14] + "_matches_mz_rt.xlsx"
    elif matching_criteria == "mz_ccs":
        output_file = input_file[:-14] + "_matches_mz_ccs.xlsx"
    elif matching_criteria == "mz_rt_ccs":
        output_file = input_file[:-14] + "_matches_mz_rt_ccs.xlsx"
    matched_df = pd.DataFrame(rows_with_matches)
    matched_df.to_excel(output_file, index=False)
    print(f"\nSuccessful analysis. View ranked matches in {output_file}.\n")

if matching_criteria in ["mz_rt_ccs_msms"]:

    print("\n...Fetching potential matches and calculating cosine similarity scores with reference MS/MS spectra...")

    # Iterate over each row in the Excel file
    for i, row in df_input.iterrows():

        # Extract the values for exact_mz, rt, ccs, ccs_calibrant, gradient, and column_type
        file_name = row["file_name"]
        sample_type = row["sample_type"]
        mz = row["mz"]
        ccs_calibrant = row["ccs_calibrant"]
        gradient = row["gradient"]
        column_type = row["column_type"]
        dt = row["dt"]
        ccs = row["ccs"]
        rt = row["rt"]
        peak_area = row["peak_area"]

        potential_matches_column = []

        # Print the current peak being analyzed to the terminal
        print(f"\n\traw file: {file_name} m/z: {mz} rt: {rt}")

        # User decides to query on m/z, rt, ccs, and MS/MS similarity score
        # Query qacs.db for matches based on selected matching criteria
            # Default criteria are: +/- 0025 Da for m/z, +/- 0.2s for rt, +/- 3% for ccs
            # Do we want to implement a m/z tolerance window?
            # The match must also have identical experimental conditions (i.e. gradient, column type, CCS calibrant)
        if matching_criteria == "mz_rt_ccs_msms":
            c.execute('''SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ?
                         AND rt BETWEEN ? AND ? AND average_ccs BETWEEN ? AND ?
                         AND gradient = ? AND ccs_calibrant = ? AND column_type = ?''', 
                      (mz - 0.025, mz + 0.025, rt - 0.2, rt + 0.2, ccs * 0.97, ccs * 1.03,
                       gradient, ccs_calibrant, column_type))

        # Get the results of the query
        matches = c.fetchall()

        # Extract potential matches from database query
        if matches:
            potential_matches = [match[0] for match in matches]
            ranked_matches = []

            # Extract MSMS spectrum for the experimental peak with potential matches
            rdr = MassLynxReader(file_name)
            if not pd.isna(dt): 
                dt = float(dt)
                m, i = rdr.get_spectrum(1, rt - 0.01, rt + 0.01, dt_min = dt - 0.1, dt_max = dt + 0.1)  # Adjust the MS function and rt, dt windows
            else:
                m, i = rdr.get_spectrum(1, rt - 0.01, rt + 0.01)

            # Convert the extracted MSMS spectrum to numpy arrays
            extracted_mz = np.array(m)
            extracted_intensity = np.array(i)

            # Filter the extracted spectrum to remove ions greater than 2 Da above the molecular ion
            extracted_peaks = list(zip(extracted_mz, extracted_intensity))
            filtered_extracted_peaks = filter_ions_by_mz(extracted_peaks, mz)

            # Normalize the extracted intensity values to the most intense peak 
            max_intensity = max(intensity for mz, intensity in filtered_extracted_peaks)
            normalized_intensity = [(mz, (intensity / max_intensity) * 100) for mz, intensity in filtered_extracted_peaks]

            # Define a threshold percentage for intensity filtering
            threshold = 0

            # Filter the extracted spectrum by intensity
            filtered_normalized_peaks = filter_peaks_by_intensity(normalized_intensity, max_intensity, threshold)

            # Split the filtered peaks back into separate arrays for further processing
            filtered_mz, filtered_intensity = zip(*filtered_normalized_peaks) if filtered_normalized_peaks else (np.array([]), np.array([]))

            similarity_list = []

            for potential_match in potential_matches:
                # User has selected qacs_experimental_msms reference database
                if database_type == "qacs_experimental_msms":
                    # Extract reference MSMS spectrum for the current potential match from qacs_experimental_msms database
                    c_msms.execute('''SELECT mz, normalized_intensity FROM experimental_msms_data 
                                          WHERE compound_name = ?''', (potential_match,))
                    msms_data = c_msms.fetchall()

                    # Deal with compounds for which there are no reference MSMS spectra in qacs_experimental_msms
                    if len(msms_data) == 0:
                        continue

                    # Pair the reference m/z and intensity values
                    reference_peaks = msms_data

                    # Calculate cosine similarity score between extracted and reference spectra
                    extracted_peaks = filtered_extracted_peaks
                    similarity = match_mz(reference_peaks, extracted_peaks, mz_tolerance, group_tolerance)
                    if np.isnan(similarity): 
                        similarity = 0
                    similarity_list.append(similarity)
 
                    # Generate a mirror plot of extracted vs. reference (experimental) MSMS spectra
                    title_mirror = "Experimental vs. Reference MS/MS Spectra \nPotential Match: {} (Score: {}) ".format(potential_match, "0" if np.isnan(similarity) else "{:.2f}".format(similarity))
                    fname_mirror = "{}/{}_{}_{}_{:.2f}_Experimental_MSMS.png".format(spectra_directory, file_name, mz, potential_match, rt)
                    mirror_plot(list(filtered_mz), list(filtered_intensity), reference_peaks, title_mirror, fname_mirror, mz_tolerance, group_tolerance)

                # User has selected qacs_theoretical_msms reference database
                elif database_type == "qacs_theoretical_msms":
                    # Extract reference MSMS spectrum for the current potential match from qacs_theoretical_msms database
                    c_theoretical_msms.execute('''SELECT mz, normalized_intensity FROM theoretical_msms_data 
                                          WHERE compound_name = ?''', (potential_match,))
                    msms_data = c_theoretical_msms.fetchall()
                    
                    # Deal with compounds for which there are no reference MSMS spectra in qacs_experimental_msms
                    if len(msms_data) == 0:
                        continue

                    # Pair the reference m/z and intensity values
                    reference_peaks = msms_data

                    # Convert the reference MSMS data to numpy arrays
                    reference_mz, reference_intensity = zip(*msms_data)
                    reference_mz = np.array(reference_mz)
                    reference_intensity = np.array(reference_intensity)

                    # Calculate cosine similarity score between extracted MSMS spectrum and reference MSMS spectrum
                    reference_mz, reference_intensity = zip(*msms_data)
                    extracted_mz, extracted_intensity = zip(*filtered_extracted_peaks)
                    reference_vector, extracted_vector = match_mz(reference_mz, reference_intensity, extracted_mz, extracted_intensity, mz_tolerance, group_tolerance)
                    similarity = similarity_score(reference_vector, extracted_vector)
                    if np.isnan(similarity): similarity = 0
                    similarity_list.append(similarity)

                    # Generate mirror plot of experimental and reference (theoretical) MSMS spectra
                    title_mirror = "Experimental vs. Reference MS/MS Spectra \nPotential Match: {} (Score: {}) ".format(potential_match, "0" if np.isnan(similarity) else "{:.2f}".format(similarity))
                    fname_mirror = "{}/{}_{}_{}_{:.2f}_Theoretical_MSMS.png".format(spectra_directory, file_name, mz, potential_match, rt)
                    mirror_plot(list(filtered_mz), list(filtered_intensity), reference_peaks, title_mirror, fname_mirror, mz_tolerance, group_tolerance)
                
            tuples = sorted(zip(similarity_list, potential_matches), reverse=True)
            ordered_potential_matches = [f"{t[1]} ({t[0]})" for t in tuples]
            ranked_matches = ', '.join(ordered_potential_matches) if ordered_potential_matches else ''
            potential_matches_column.append(", ".join(potential_matches))
            row_data = [file_name, sample_type, mz, ccs_calibrant, gradient, column_type, dt, ccs, rt, peak_area, potential_matches_column, ranked_matches]
            final_rows.append(row_data)

    # Create a new dataframe from the final rows
    columns = ["file_name", "sample_type", "mz", "ccs_calibrant", "gradient", "column_type", "dt", "ccs", "rt", "peak_area", "potential_matches", "ranked_matches"]
    df_output = pd.DataFrame(final_rows, columns=columns)

    # Write the dataframe to an Excel file
    if database_type == "qacs_experimental_msms":
        output_file = input_file[:-14] + "_matches_ranked_experimental.xlsx"
    elif database_type == "qacs_theoretical_msms":
        output_file = input_file[:-14] + "_matches_ranked_theoretical.xlsx"
    
    # Handle the ranked_matches column 
    def format_ranked_matches(value):
        try:
            tuples = value.split(", ")  # Split the potential matches and scores
            matches = [t.split(" (")[0] for t in tuples]  # Extract the potential matches
            scores = [t.split(" (")[1].split(")")[0] for t in tuples]  # Extract the scores
            scores = ["0" if score == "nan" else score for score in scores]
            rounded_scores = [str(round(float(score), 2)) for score in scores]  # Round the scores to two decimal places
            formatted_matches = [f"{match} ({score})" for match, score in zip(matches, rounded_scores)]  # Combine matches and scores
            return ", ".join(formatted_matches)
        except (ValueError, IndexError):
            return value
            
    # Handle the potential_matches column
    def format_potential_matches(value):
        if isinstance(value, list):
            unique_matches = list(set(value)) # Remove duplicate potential matches
            return ", ".join(unique_matches)
        return value

    # Export the dataframe to Excel
    df_output["ranked_matches"] = df_output["ranked_matches"].apply(format_ranked_matches)  # Display 2 decimal places in ranked_matches column
    df_output["potential_matches"] = df_output["potential_matches"].apply(format_potential_matches).str.replace(",", ", ")  # Add space after comma separator
    df_output.to_excel(output_file, index=False)
    print(f"\nSuccessful analysis. View ranked matches in {output_file}.\n")

# Close the database connections
conn.close()
conn_msms.close()
conn_theoretical_msms.close()