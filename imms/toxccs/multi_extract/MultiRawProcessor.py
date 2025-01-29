"""
toxccs/multi_extract/MultiRawProcessor.py

Ryan Nguyen (ryan97@uw.edu)
12/21/24

description:
        Main logic for extracting and processing Waters .raw LC-IM-MS/MS data. Specific logic included for handling complex ion mobilogram profiles.
"""

import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from dhrmasslynxapi.reader import MassLynxReader
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRawXL
from toxccs.utils.gaussian import process_multi_chromatogram, gaussian_smooth_pick
from toxccs.utils.peaks import observed_mz
from toxccs.utils.plotting import combined_multi_figure
from toxccs.utils.metrics import calculate_peak_resolution


# Set global font conditions for figures
params = {"font.size": 8, "font.family": "Arial", "font.weight": "bold"}
plt.rcParams.update(params)


class MultiRawProcessor:
    def __init__(self, target_list, reference_is_file=None):
        """
        MultiRawProcessor.__init__
        description:
                Initializes a new MultiRawProcessor with a target list containing precursor ions and paths to Waters .raw data files.
        parameters:
                target_list (str) -- path to the Excel file (.xlsx) containing the target list.
                reference_is_list (str) -- path to the Excel file (.xlsx) containing the reference internal standard list (OPTIONAL).
        """

        self.target_list = target_list
        self.reference_is_df = None
        if reference_is_file:
            self.read_reference_is(reference_is_file)
        self.output_rows = []

    def export_to_excel(self):
        """
        MultiRawProcessor.export_to_excel
        description:
                Exports the extracted data to an Excel file (.xlsx). Specific data included for complex mobilogram profiles.
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
                "LC FWHM Flag",
                "EIM FWHM (ms)",
                "Two-Peak EIM Resolution",
            ],
        )

        # Name the output file
        base_name = os.path.splitext(os.path.basename(self.target_list))[0]
        output_file = "{}_processed.xlsx".format(base_name)
        df_output.to_excel(output_file, index=False)
        print("\nProcessing successful. View results in {}.".format(output_file))

        return output_file

    def read_reference_is(self, reference_is_file):
        """
        MultiRawProcessor.read_reference_is
        description:
                Reads and stores the reference internal standards Excel file (.xlsx). Intended to be used primarily for retention time alignment and/or quantitation of target compounds.
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

    def multi_extract(
        self, calibration_file, ms1_function, mobility_function, mz_tolerance
    ):
        """
        MultiRawProcessor.multi_extract
        description:
                Main execution sequence for extracting and processing Waters .raw LC-IM-MS/MS data. Specific logic for handling complex ion mobilogram profiles.
        parameters:
                calibration_file (str) -- path to the CCS calibration file.
                ms1_function (int) -- MS1 function number.
                mobility_function (int) -- mobility function number.
                mz_tolerance (float) -- Da tolerance around the target m/z value to consider. Default is 0.025.
        returns:
                (.xlsx) -- Excel spreadsheet containing the extracted data for each feature.
                (.png) -- Image displaying the chemical structure, MS1 scan, extracted ion chromatogram (EIC), and extracted ion mobilogram (EIM).
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

        # Generate a new folder called "Extracted Data" in the current directory if not already created
        data_directory = "Extracted Reprocessed Data"
        if not os.path.exists(data_directory):
            os.makedirs(data_directory)

        # Begin main execution sequence for extraction
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

            # Begin sequence for identifying monoisotopic peak within the MS1 can
            # Average all scans within the acquisition time
            # Default is 0.2 to 1.15 min
            rdr = MassLynxReader(file_name)
            mz_spectrum, intensity_spectrum = rdr.get_spectrum(ms1_function, 0.2, 1.15)

            # Smooth and fit the extracted MS1 spectrum using Gaussian convolution
            # This process is analogous to centroiding the profile data, as outlined in MSnbase
            identified_peaks, smoothed_intensity = gaussian_smooth_pick(
                mz_spectrum, intensity_spectrum
            )

            # Identify the monoisotopic peak in the MS1 centroided data
            (monoisotopic_mz, monoisotopic_intensity), _ = observed_mz(
                identified_peaks, mz, self.mz_tolerance
            )
            if monoisotopic_mz == 0:
                continue

            # Begin sequence for extracting ion chromatogram
            # Extract m/z-selected ion chromatogram (EIC) from the .raw file using identified monoisotopic m/z value
            # Requires the correct MS1 function number and desired m/z tolerance
            rdr = MassLynxReader(file_name)
            rt, rt_i = rdr.get_chrom(
                self.ms1_function, float(monoisotopic_mz), self.mz_tolerance
            )

            # Store extracted rt, rt_i data as arrays for further processing
            rt = np.array(rt)
            rt_i = np.array(rt_i)

            # Smooth and fit EIC with multi-Gaussian function
            (
                valid_peak_indices,
                rt_list,
                eic_areas,
                peak_ranges,
                _,
                mu_values,
                fwhm_flags,
                rt_fwhm_values,
            ) = process_multi_chromatogram(rt, rt_i, apply_fifty_percent_rule=False)

            # Begin sequence for extracting ion mobilogram
            # Iterate over each extracted LC ion chromatogram peak
            for (
                rt_value,
                eic_area,
                (start_idx, end_idx),
                mu,
                fwhm_flag,
                rt_fwhm_value,
            ) in zip(
                rt_list, eic_areas, peak_ranges, mu_values, fwhm_flags, rt_fwhm_values
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
                    print(
                        f"Adjusted RT range due to close spacing: {rt_start} to {rt_end}"
                    )

                # Extract (m/z,rt)-selected ion mobilogram (EIM) for each identified LC peak
                t, dt_i = rdr.get_filtered_chrom(
                    self.mobility_function,
                    float(monoisotopic_mz),
                    self.mz_tolerance,
                    rt_min=rt_start,
                    rt_max=rt_end,
                )
                dt = None

                # Store extracted t, dt_i data as arrays for further processing
                t = np.array(t)
                dt_i = np.array(dt_i)

                # Smooth and fit EIM with multi-Gaussian function
                (
                    _,
                    _,
                    _,
                    _,
                    _,
                    eim_mu_values,
                    _,
                    eim_fwhm_values,
                ) = process_multi_chromatogram(t, dt_i, apply_fifty_percent_rule=True)

                # Convert drift times to CCS values and store results
                ccs_values = [
                    round(self.cal_data.calibrated_ccs(monoisotopic_mz, eim_mu), 2)
                    for eim_mu in eim_mu_values
                ]

                # Calculate the EIC peak height (i.e., max intensity)
                peak_height = max(rt_i[start_idx : end_idx + 1])

                # Append the extracted data to output_rows list
                mobilogram_data = []
                for eim_mu, ccs, eim_fwhm in zip(
                    eim_mu_values, ccs_values, eim_fwhm_values
                ):
                    mobilogram_data.append(
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
                            eim_mu,
                            ccs,
                            rt_value,
                            peak_height,
                            eic_area,
                            fwhm_flag,
                            eim_fwhm,
                            None,
                        ]
                    )

                # Calculate the resolution for multiple EIM peaks at the same retention time
                # This calculation could be moved to a helper function
                if len(mobilogram_data) > 1:
                    for i in range(len(mobilogram_data) - 1):
                        dt_a = mobilogram_data[i][10]
                        fwhm_a = mobilogram_data[i][16]
                        dt_b = mobilogram_data[i + 1][10]
                        fwhm_b = mobilogram_data[i + 1][16]

                        # Skip resolution calculation if any of the values are missing or invalid
                        if any(
                            v in [None, "", "NaN"] for v in [dt_a, fwhm_a, dt_b, fwhm_b]
                        ):
                            continue

                        # Calculate resolution
                        resolution = calculate_peak_resolution(
                            dt_a, fwhm_a, dt_b, fwhm_b
                        )

                        # Store the resolution in the row corresponding to the lower drift time
                        if resolution is not None:
                            if dt_a < dt_b:
                                mobilogram_data[i][-1] = resolution
                            else:
                                mobilogram_data[i + 1][-1] = resolution
                self.output_rows.extend(mobilogram_data)

                # Generate file name for combined figure
                fname_combined = f"{data_directory}/{compound_name}_{mz}_{adduct}_{rt_value}_{file_name}.png"

                # Generate combined figure
                combined_multi_figure(
                    rt,
                    rt_i,
                    t,
                    dt_i,
                    file_name,
                    compound_name,
                    adduct,
                    mz,
                    monoisotopic_mz,
                    monoisotopic_intensity,
                    rt_value,
                    mu_values,
                    fname_combined,
                    mz_spectrum,
                    intensity_spectrum,
                    smiles,
                )

        # Export DataFrame containing extracted spectral features to Excel file (.xlsx)
        self.output_file = self.export_to_excel()
