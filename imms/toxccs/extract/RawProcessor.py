"""
toxccs/extract/RawProcessor.py

Ryan Nguyen (ryan97@uw.edu)
12/20/24

description:
        Main logic for extracting and processing Waters .raw LC-IM-MS/MS data. 
"""

import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from matplotlib import rcParams
from dhrmasslynxapi.reader import MassLynxReader
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRawXL
from toxccs.utils.gaussian import (
    process_chromatogram,
    peak_fit,
    gaussian_smooth_pick,
    gaussian_fit,
    peak_fit,
)
from toxccs.utils.peaks import observed_mz
from toxccs.utils.plotting import combined_figure
from toxccs.utils.threshold import fwhm_threshold


# Set global font conditions for figures
params = {"font.size": 8, "font.family": "Arial", "font.weight": "bold"}
plt.rcParams.update(params)


class RawProcessor:
    def __init__(self, target_list, reference_is_file=None):
        """
        RawProcessor.__init__
        description:
                Initializes a new RawProcessor with a target list containing precursor ions and paths to Waters .raw data files.
        parameters:
                target_list (str) -- path to the Excel file (.xlsx) containing the target list.
                reference_is_list (str) -- path to the Excel file (.xlsx) containing the reference internal standard list (OPTIONAL).
        """

        self.target_list = target_list
        self.reference_is_df = None
        if reference_is_file:
            self.read_reference_is(referene_is_file)
        self.output_rows = []

    def export_to_excel(self):
        """
        RawProcessor.export_to_excel
        description:
                Exports the extracted data to an Excel file (.xlsx).
        returns:
                (str) -- path to the output Excel file (xlsx).
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

        # Name the output file
        base_name = os.path.splitext(os.path.basename(self.target_list))[0]
        output_file = "{}_processed.xlsx".format(base_name)
        df_output.to_excel(output_file, index=False)

        # Print statement upon successful extraction sequence
        print("\nProcessing successful. View results in {}.".format(output_file))

        return output_file

    def read_reference_is(self, reference_is_file):
        """
        RawProcessor.read_reference_is
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

    def extract(self, calibration_file, ms1_function, mobility_function, mz_tolerance):
        """
        RawProcessor.extract
        description:
                Main execution sequence for extracting and processing Waters .raw LC-IM-MS/MS data.
        parameters:
                calibration_file (str) -- path to the CCS calibration file.
                ms1_function (int) -- MS1 function number.
                mobility_function (int) -- mobility function number.
                mz_tolerance (float) -- Da tolerance around the target m/z value to consider. Default is 0.025 Da.
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
        data_directory = "Extracted Data"
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

            # Extract the path to .raw file, target m/z value, sample type (i.e., sample or blank), gradient, adduct type, compound name, and SMILES string
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

            # Begin sequence for identifying monoisotopic peak within the MS1 scan
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

            # Begin sequence for extracting ion mobilogram
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

                # Smooth and fit EIM to a single Gaussian function
                A, B, C = peak_fit(t, dt_i)
                t_refined = np.arange(min(t), max(t), 0.01)
                fit_i = gaussian_fit(t_refined, A, B, C)

                # Apply full width at half maximum (FWHM) and EIM intensity thresholds
                if fwhm_threshold(C, dt_i):
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
                # This calculation could be moved to a helper function
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

                # Calculate FWHM for output figure
                # This calculation could also be pulled from helper function
                fwhm = C * 2.355

                # Generate file name for combined figure
                fname_combined = f"{data_directory}/{compound_name}_{mz}_{adduct}_{rt_value}_{file_name}.png"

                # Generate combined figure
                combined_figure(
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
