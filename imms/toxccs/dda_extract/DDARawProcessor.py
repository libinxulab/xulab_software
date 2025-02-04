"""
toxccs/dda_extract/DDARawProcessor.py

Ryan Nguyen (ryan97@uw.edu)
2/4/25

description:
        Main logic for extracting MS/MS fragmentation spectra from Waters Fast DDA .raw files. Note that the data must be centroided prior to implementation.
"""

import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from dhrmasslynxapi.reader import MassLynxReader
from toxccs.utils.gaussian import process_chromatogram
from toxccs.utils.plotting import combined_dda_figure
from toxccs.utils.peaks import observed_mz
from toxccs.utils.metrics import calculate_mass_error


# Set global font conditions for figures
params = {"font.size": 8, "font.family": "Arial", "font.weight": "bold"}
plt.rcParams.update(params)


class DDARawProcessor:
    def __init__(self, target_list):
        """
        DDARawProcessor.__init__
        description:
                Initializes a new RawProcessor with a target list containing precursor ions and paths to Waters .raw data files.
        parameters:
                target_list (str) -- path to the Excel file (.xlsx) containing the target list.
        """

        self.target_list = target_list
        self.reference_is_df = None
        self.output_rows = []

    def export_to_excel(self):
        """
        DDARawProcessor.export_to_excel
        description:
                Exports the extracted MS/MS fragmentation spectrum to an Excel file (.xlsx).
        returns:
                (str) -- path to the output Excel file (xlsx).
        """

        for df_output, output_name in self.output_rows:
            df_output.to_excel(
                output_name,
                index=False,
                columns=["m/z", "Intensity"],
                header=["m/z", "Intensity"],
            )

    def export_features_to_excel(self, feature_data, output_file):
        """
        DDARawProcessor.export_features_to_excel
        description:
                Exports the extracted spectral features (i.e., observed m/z and retention time) to an Excel file (.xlsx).
        parameters:
                feature_data (list) -- list of extracted and processed spectral features.
                output_file (str) -- path to the output Excel file (.xlsx).
        """
        df_features = pd.DataFrame(
            feature_data,
            columns=[
                "File Name",
                "Compound Name",
                "SMILES",
                "Adduct",
                "Target m/z",
                "Observed m/z",
                "Channel Number",
                "Mass Error (ppm)",
                "Sample Type",
                "Gradient",
                "Column Type",
                "Observed Retention Time (min)",
                "Flags",
            ],
        )
        df_features.to_excel(output_file, index=False)
        print(f"\nProcessing successful. View results in {output_file}.")

    def dda_extract(self, mz_tolerance, num_precursors):
        """
        DDARawProcessor.dda_extract
        description:
                Main sequence for extracting MS/MS fragmentation spectra from Waters DDA .raw files.
        parameters:
                mz_tolerance (float) -- Da tolerance around the target m/z value to consider. Default is 0.025 Da.
                num_precursors (int) -- Number of precursor ions subjected to fragmentation at each MS/MS event.
        returns:
                (.xlsx) -- Excel spreadsheet containing the extracted MS/MS fragmentation spectrum (i.e., list of paired m/z and intensity values) for each input target.
                (.xlsx) -- Single Excel spreadsheet containing the extracted features for all targets in the input precursor list.
                (.png) -- Image displaying the chemical structure, MS1 scan, extracted ion chromatogram (EIC), and extracted MS/MS fragmentation spectrum.
        """

        # Set the input parameters
        self.mz_tolerance = mz_tolerance
        self.num_precursors = num_precursors

        # Initialize output rows for exported data
        feature_data = []
        self.output_rows = []

        # Read target list into a pandas DataFrame
        df_input = pd.read_excel(self.target_list)
        df_input.reset_index(drop=True, inplace=True)

        # Generate a new folder called "Extracted Spectra" in the current directory if not already created
        data_directory = "Extracted Fragmentation Spectra"
        if not os.path.exists(data_directory):
            os.makedirs(data_directory)

        # Begin main sequence for extraction
        # Iterate over each row in the target list DataFrame
        print(
            "\n...Extracting MS/MS fragmentation spectra from the raw files using precursor target list {}...".format(
                self.target_list
            )
        )
        for i, row in df_input.iterrows():

            # Extract the path to .raw file, target m/z value, sample type (i.e., sample or blank), gradient, adduct type, compound name, and SMILES string
            (
                file_name,
                mz,
                sample_type,
                gradient,
                adduct,
                compound_name,
                smiles,
                column_type,
            ) = (
                row["File Name"],
                row["Target m/z"],
                row["Sample Type"],
                row["Gradient"],
                row["Adduct"],
                row["Compound Name"],
                row["SMILES"],
                row["Column Type"],
            )

            # Print the current m/z being queried to the terminal
            print(
                "\n\traw file: {} compound: {} m/z: {}".format(
                    file_name, compound_name, mz
                )
            )

            # Check if the file exists before proceeding
            if not os.path.exists(file_name):
                print(f"\tFile not found: {file_name}.")
                continue

            # Begin sequence for extracting MS1-level ion chromatogram
            # Extract m/z-selected ion chromatogram (EIC) from the .raw file using theoretical m/z value of target and specified tolerance
            # MS1 survey scan data are always stored in Channel 1
            rdr = MassLynxReader(file_name)
            rt, rt_i = rdr.get_chrom(0, float(mz), self.mz_tolerance)

            # Store extracted rt, rt_i data as arrays for further processing
            rt = np.array(rt)
            rt_i = np.array(rt_i)

            # Smooth and fit EIC with multi-Gaussian function
            (
                peak_indices,
                rt_list,
                areas,
                peak_ranges,
                popt_multi_gaussian,
                mu_values,
            ) = process_chromatogram(
                rt,
                rt_i,
                file_name,
                compound_name,
                adduct,
                mz,
                sample_type,
                generate_images=True,
            )

            # Dynamically calculate EIC time bounds
            found_peak = False

            # Identify most intense EIC if multiple peaks are present in MS1 survey scan for a given target precursor
            if len(rt_list) > 0:
                max_intensity_index = np.argmax([rt_i[peak] for peak in peak_indices])
                peak_apex = rt_list[max_intensity_index]

                # Use peak indices to define the window
                start_idx, end_idx = peak_ranges[max_intensity_index]
                rt_min = rt[start_idx]
                rt_max = rt[end_idx]

                # Handle cases where peak indices are determined to be equal due to close spacing of raw data
                # Default to use a fixed 0.3 min window (i.e., +/- 0.15 min around the apex)
                if rt_min == rt_max:
                    fixed_bound = 0.15
                    rt_min = max(rt[0], peak_apex - fixed_bound)
                    rt_max = min(rt[-1], peak_apex + fixed_bound)
                found_peak = True
            else:
                print("\t\tNo valid MS1 chromatographic peak identified for target.")
                continue

            # Extract MS1 survey scan scan coordinates defined by rt_min and rt_max
            try:
                scan_indices_ms1 = rdr._MassLynxReader__get_scan_indices(
                    0, rt_min, rt_max
                )
                if scan_indices_ms1:
                    scan_min, scan_max = scan_indices_ms1[0], scan_indices_ms1[-1]
                    scan_times = [rdr.scan_times[0][idx] for idx in scan_indices_ms1]
                    # print(f"\tMS1 Scan Range: {scan_min} - {scan_max}")
                    # print("\tScan Index - Retention Time Mapping:")
                    # for scan_idx, scan_time in zip(scan_indices_ms1, scan_times):
                    # print(f"\tScan {scan_idx}: {scan_time:.4f} min")
                else:
                    print(
                        "\t\tNo corresponding MS1 scans found within the specified retention time window."
                    )
            except Exception as e:
                print(f"\t\tError retrieving MS1 scan numbers: {e}")

            # Initialize values and lists to track the precursor ion peak across MS1 survey scans
            most_intense_ms1_scan_index = None
            most_intense_precursor_mz = None
            most_intense_precursor_intensity = -np.inf
            best_top_peaks_mz = []
            best_top_peaks_intensity = []

            # Extract the MS1 survey scan spectrum for each scan in scan_indices_ms1
            try:
                for scan_idx in scan_indices_ms1:
                    scan_mz, scan_intensity = rdr.scan_reader.ReadScan(0, scan_idx)

                    # Save the extracted MS1 survey spectrum from each scan
                    scan_mz = np.array(scan_mz)
                    scan_intensity = np.array(scan_intensity)

                    # Identify the precursor ion peak within the extracted MS1 survey spectrum
                    (
                        observed_precursor_ms1_mz,
                        observed_precursor_ms1_intensity,
                    ), _ = observed_mz(
                        zip(
                            scan_mz,
                            scan_intensity,
                        ),
                        mz,
                        self.mz_tolerance,
                    )

                    # Identify the top precursor_num most intense peaks in the scan
                    top_indices = np.argsort(scan_intensity)[-num_precursors:][::-1]
                    top_peaks_mz = scan_mz[top_indices]
                    top_peaks_intensity = scan_intensity[top_indices]

                    # Track the scan with the most intense precursor peak
                    if (
                        observed_precursor_ms1_intensity
                        > most_intense_precursor_intensity
                    ):
                        most_intense_precursor_intensity = (
                            observed_precursor_ms1_intensity
                        )
                        most_intense_precursor_mz = observed_precursor_ms1_mz
                        most_intense_ms1_scan_index = scan_idx
                        best_top_peaks_mz = top_peaks_mz
                        best_top_peaks_intensity = top_peaks_intensity

                # Print extracted data
                if most_intense_ms1_scan_index is not None:
                    print(
                        f"\n\t\tMost intense MS1 survey scan: {most_intense_ms1_scan_index} ({rdr.scan_times[0][most_intense_ms1_scan_index]:.4f} min)"
                    )
                    print(
                        f"\t\tPrecursor Ion Peak: m/z {most_intense_precursor_mz:.4f} Intensity {most_intense_precursor_intensity:.1f}"
                    )
                    print(f"\t\tTop {num_precursors} Peaks in Most Intense MS1 Scan:")
                    print("")
                    for rank, (mz_val, intensity_val) in enumerate(
                        zip(best_top_peaks_mz, best_top_peaks_intensity), start=1
                    ):
                        print(
                            f"\t\t\t{rank}. m/z: {mz_val:.4f} Intensity: {intensity_val:.1f}"
                        )
                    print("")
                else:
                    print("\n\t\tNo valid MS1 scan with intense precursor ion found.")
            except Exception as e:
                print(f"\n\t\tError processing MS1 survey scans: {e}")

            # Begin sequence for extracting MS/MS fragmentation spectrum
            # In a typical Waters Fast-DDA .raw file, fragmentation data from the most intense precursor at a given MS/MS timepoint are store in Channel 2 (MS1 survey scan data are stored in Channel 1)
            # Data from next most intense precursor are stored in Channel 3, followed by Channel 4, etc.
            # The last Channel contains LockMass data, if included
            # Identify and extract the MS/MS fragmentation spectrum with the most intense target precursor ion
            max_intensity = -np.inf
            correct_channel = None
            ms2, ms2_i = None, None
            nonzero_intensity_channels = []

            # Iterate through the possible Channels
            for channel_num in range(1, num_precursors + 1):
                try:

                    # Extract MS/MS fragmentation spectrum from the current Channel
                    # rt_m = float(peak_apex) + 0.01 optimized during testing
                    spectrum_mz, spectrum_intensity = rdr.get_spectrum(
                        channel_num, rt_min, rt_max
                    )
                    spectrum_mz = np.array(spectrum_mz)
                    spectrum_intensity = np.array(spectrum_intensity)

                    # Identify the precursor ion peak within the MS/MS fragmentation spectrum
                    (observed_precursor_mz, observed_precursor_intensity), _ = (
                        observed_mz(
                            zip(spectrum_mz, spectrum_intensity), mz, self.mz_tolerance
                        )
                    )

                    # Save the Channel number for all instances where the precursor ion peak is detected
                    if observed_precursor_intensity > 0:
                        nonzero_intensity_channels.append(channel_num)

                    # Update the correct Channel if this spectrum has a higher precursor ion intensity
                    # The assumption is that the Channel storing the desired fragmentation data will have the most intense target precursor peak
                    # All other channels should have target precursor intensities close to 0
                    if observed_precursor_intensity > max_intensity:
                        max_intensity = observed_precursor_intensity
                        correct_channel = channel_num

                        # Save the correct observed monoisotopic peak for mass error calculation
                        ms2_monoisotopic_mz = observed_precursor_mz

                        # Save the correct function number adjusted for indexing
                        adj_correct_channel = correct_channel + 1

                except Exception as e:
                    print(f"\n\t\tChannel {channel_num} could not be processed: {e}")
            if nonzero_intensity_channels:
                formatted_channels = ", ".join(
                    map(str, [ch + 1 for ch in nonzero_intensity_channels])
                )
                print(f"\t\tCandidate MS/MS Channels Identified: {formatted_channels}")
            if correct_channel is None:
                print("\n\t\tNo valid MS2 Channel found for target.")
                continue

            # Iterate through individual scans in the accumulation window to find the scan containing the most intense target ion precursor
            scan_indices = rdr._MassLynxReader__get_scan_indices(
                correct_channel, rt_min, rt_max
            )
            most_intense_scan = None
            most_intense_value = -np.inf
            for scan_index in scan_indices:
                scan_mz, scan_intensity = rdr.scan_reader.ReadScan(
                    correct_channel, scan_index
                )
                scan_mz = np.array(scan_mz)
                scan_intensity = np.array(scan_intensity)

                # Identify the precursor ion peak within the extracted MS/MS fragmentation spectrum
                (observed_precursor_mz, observed_precursor_intensity), _ = observed_mz(
                    zip(scan_mz, scan_intensity), mz, self.mz_tolerance
                )
                if observed_precursor_intensity > most_intense_value:
                    most_intense_value = observed_precursor_intensity
                    most_intense_scan = (scan_mz, scan_intensity)
                    most_intense_scan_index = scan_index
                if most_intense_scan is not None:
                    ms2, ms2_i = most_intense_scan

                # Implement MS1 survey scan intensity check to validate correct_channel
                # Retrieve MS1 spectrum from scan containing the most intense precursor ion peak determined on line 276
                # Check if the precuror ion is among the top precursor_num most intense peaks within the MS1 survey scan
                flag_status = ""
                needs_precursor_check = (
                    most_intense_precursor_mz not in best_top_peaks_mz
                )

                # Define updated variables for print statements and data export
                final_channel = None
                final_intensity = None
                final_scan_index = None
                observed_precursor_mz_new = None

                if needs_precursor_check:

                    # If there is only one valid channel_num in nonzero_intensity_channels, flag the compound in output
                    if len(nonzero_intensity_channels) == 1:
                        flag_status = "FLAG"
                        final_channel = nonzero_intensity_channels[0]

                    # If there are exactly two valid channel_num in nonzero_intensity_channels, switch to the second Channel and apply scan iteration logic
                    elif len(nonzero_intensity_channels) == 2:
                        new_channel = nonzero_intensity_channels[1]
                        correct_channel = new_channel

                        # Apply scan iteration logic within new correct_channel
                        scan_indices_new = rdr._MassLynxReader__get_scan_indices(
                            new_channel, rt_min, rt_max
                        )
                        most_intense_value_new = -np.inf
                        for scan_index in scan_indices_new:
                            scan_mz, scan_intensity = rdr.scan_reader.ReadScan(
                                new_channel, scan_index
                            )
                            scan_mz = np.array(scan_mz)
                            scan_intensity = np.array(scan_intensity)

                            # Identify the precursor ion peak with the extracted MS/MS fragmentation spectrum
                            (
                                observed_mz_temp,
                                observed_intensity_temp,
                            ), _ = observed_mz(
                                zip(scan_mz, scan_intensity), mz, self.mz_tolerance
                            )

                            # Update the correct scan if this spectrum has a higher precursor ion intensity
                            if observed_intensity_temp > most_intense_value_new:
                                most_intense_value_new = observed_intensity_temp
                                final_scan_index = scan_index
                                observed_precursor_mz_new = observed_mz_temp
                        final_channel = new_channel
                        final_intensity = most_intense_value_new

                    # If there are more than 3 valid channel_num in nonzero_intensity_channels, evaluate these
                    elif len(nonzero_intensity_channels) >= 3:
                        best_channel = None
                        best_precursor_intensity = -np.inf
                        for new_channel in nonzero_intensity_channels[1:3]:
                            scan_indices_new = rdr._MassLynxReader__get_scan_indices(
                                new_channel, rt_min, rt_max
                            )
                            most_intense_value_new = -np.inf
                            most_intense_scan_index_new = None
                            for scan_index in scan_indices_new:

                                # Apply scan iteration logic within MS/MS fragmentation Channels
                                scan_mz, scan_intensity = rdr.scan_reader.ReadScan(
                                    new_channel, scan_index
                                )
                                scan_mz = np.array(scan_mz)
                                scan_intensity = np.array(scan_intensity)

                                # Identify the precursor ion peak within the extracted MS/MS fragmentation spectrum
                                (
                                    observed_mz_temp,
                                    observed_intensity_temp,
                                ), _ = observed_mz(
                                    zip(scan_mz, scan_intensity), mz, self.mz_tolerance
                                )

                                # Update the correct scan if this spectrum has a higher precursor ion intensity
                                if observed_intensity_temp > most_intense_value_new:
                                    most_intense_value_new = observed_intensity_temp
                                    most_intense_scan_index_new = scan_index
                                    temp_mz = observed_mz_temp

                            # Choose the better Channel based on precuror intensty
                            if most_intense_value_new > best_precursor_intensity:
                                best_precursor_intensity = most_intense_value_new
                                best_channel = new_channel
                                best_observed_mz = temp_mz
                        final_channel = best_channel
                        final_intensity = best_precursor_intensity
                        observed_precursor_mz_new = best_observed_mz
            print(observed_precursor_mz_new)
            if needs_precursor_check:
                print(
                    f"\t\t...Precursor ion peak NOT among {num_precursors} most intense peaks in MS1 surey scan..."
                )
                if len(nonzero_intensity_channels) == 2:
                    print(
                        f"\t\tEvaluating MS/MS fragmentation data from {nonzero_intensity_channels[1]}."
                    )
                else:
                    print(
                        f"\t\tEvaluating MS/MS fragmentation data from additional valid Channels."
                    )
            else:
                print(
                    f"\t\tPrecursor ion peak IS among {num_precursors} most intense peaks in MS1 survey scan..."
                )
                print(f"\t\tEvaluating original Channel.")

            # Generate names for output spectrum files
            base_name = f"{compound_name}_{mz}_{adduct}_{os.path.splitext(os.path.basename(file_name))[0]}"
            spectrum_output_path = os.path.join(
                data_directory, f"{base_name}.raw_MSMS.xlsx"
            )
            figure_output_path = os.path.join(data_directory, f"{base_name}.raw.png")

            # Remove extraneous peaks greater than precursor target m/z
            filtered_ms2 = ms2[ms2 < observed_precursor_mz + 4]
            filtered_ms2_i = ms2_i[ms2 < observed_precursor_mz + 4]

            # Initialize a fallback for the observed precursor ion if filtered spectrum is empty
            observed_precursor_fallback_mz = None
            observed_precursor_fallback_intensity = None

            # Handle case where filtered MS2 spectrum is empty
            if len(filtered_ms2) == 0 or len(filtered_ms2_i) == 0:
                print(
                    "\tNo MS2 data found after mass filtering. Attempting to fetch data from MS1 survey scan."
                )
                try:
                    ms1_mz, ms1_intensity = rdr.get_spectrum(0, rt_min, rt_max)
                    ms1_mz = np.array(ms1_mz)
                    ms1_intensity = np.array(ms1_intensity)
                    (
                        observed_precursor_fallback_mz,
                        observed_precurspr_fallback_intensity,
                    ), _ = observed_mz(
                        zip(ms1_mz, ms1_intensity), mz, self.mz_tolerance
                    )
                except Exception as e:
                    print(f"\tError fetching MS1 precuror ion: {e}")

            # calculate mass error in ppm using fallback value if MS2 spectrum is empty
            observed_mz_final = (
                observed_precursor_fallback_mz
                if observed_precursor_fallback_mz
                else (
                    observed_precursor_mz_new
                    if needs_precursor_check and observed_precursor_mz_new is not None
                    else ms2_monoisotopic_mz
                )
            )

            # Calculate value
            mass_error = (
                calculate_mass_error(mz, observed_mz_final)
                if observed_mz_final is not None
                else None
            )

            # Set "MS1" if the observed precursor ion peak was extracted from the MS1 survey scan
            channel_output = (
                "MS1"
                if observed_mz_final == observed_precursor_fallback_mz
                else (
                    final_channel + 1
                    if final_channel is not None
                    else adj_correct_channel
                )
            )

            # Append the extracted data to feature_data list
            feature_data.append(
                [
                    file_name,
                    compound_name,
                    smiles,
                    adduct,
                    mz,
                    observed_mz_final,
                    channel_output,
                    mass_error,
                    sample_type,
                    gradient,
                    column_type,
                    peak_apex,
                    flag_status,
                ]
            )

            # Append the extracted MS/MS fragmentation spectrum to output_rows
            self.output_rows.append(
                (
                    pd.DataFrame({"m/z": filtered_ms2, "Intensity": filtered_ms2_i}),
                    spectrum_output_path,
                )
            )

            # Generate combined figure
            combined_dda_figure(
                rt,
                rt_i,
                filtered_ms2,
                filtered_ms2_i,
                file_name,
                compound_name,
                adduct,
                mz,
                figure_output_path,
                smiles,
                rt_min,
                rt_max,
            )

        # Export feature data
        base = os.path.splitext(os.path.basename(self.target_list))[0]
        feature_output_file = f"{base}_processed.xlsx"
        self.export_features_to_excel(feature_data, feature_output_file)

        # Export fragmentation spectrum data
        self.export_to_excel()
