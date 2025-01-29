"""
    toxccs/utils/peaks.py

    Ryan Nguyen (ryan97@uw.edu)
    12/20/24

    description:
            Module with utilities to process mass spectrometric peaks. 
"""

import numpy as np


def subtract_background_spectrum(
    main_mz, main_intensity, background_mz, background_intensity, mz_tolerance
):
    """
    utils/subtract_background_spectrum
    description:
            Subtracts overlapping peaks from the extracted MS/MS spectrum based on a background spectrum.
    parameters:
            main_mz (numpy.ndarray) -- m/z values of the main MS/MS spectrum.
            main_intensity (numpy.ndarray) -- intensity values of the main MS/MS spectrum.
            background_mz (numpy.ndarray) -- m/z values of the background spectrum.
            background_intensity (numpy.ndarray) -- intensity values of the background spectrum.
            mz_tolerance (float) -- tolerance for matching peaks in Da.
    returns:
            filtered_mz (numpy.ndarray) -- filtered m/z values after background subtraction.
            filtered_intensity (numpy.ndarray) -- filtered intensity values after background subtraction.
    """
    if len(background_mz) == 0:
        return main_mz, main_intensity

    # Create a mask to exclude peaks that overlap with the background spectrum
    mask = np.array(
        [
            not any(
                abs(main_peak - bg_peak) <= mz_tolerance for bg_peak in background_mz
            )
            for main_peak in main_mz
        ]
    )

    # Filter the main spectrum
    filtered_mz = main_mz[mask]
    filtered_intensity = main_intensity[mask]

    return filtered_mz, filtered_intensity


def monoisotopic_peak(identified_peaks, target_mz, tolerance=0.025):
    """
    utils/monoisotopic_peak
    description:
            Identifies the monoisotopic peak within the target rt-selected MS1 spectrum. The current implemenation only considers the M+0 peak because all [M]+ analytes carry a permanent positive charge.
    parameters:
            identified_peaks (tuple) -- extracted and centroided MS1 spectrum.
            target_mz (float) -- target precursor m/z value to compare with extracted peaks.
            tolerance (float) -- Da tolerance around the target m/z value to consider. Default is 0.025 Da.
    returns:
            (float) -- extracted m/z value of the monoisotopic peak.
    """

    # Identify peaks from the extracted MS1 function that fall within the M+3 window
    potential_peaks = [
        (mz, intensity)
        for mz, intensity in identified_peaks
        if target_mz - tolerance <= mz <= target_mz + tolerance
    ]

    # If no peaks are found within tolerance, return the target m/z
    if not potential_peaks:

        return target_mz

    # Find the peak with the highest intensity within the tolerance window
    highest_intensity_peak = max(potential_peaks, key=lambda x: x[1])

    return highest_intensity_peak[0]


def observed_mz(identified_peaks, target_mz, mz_tolerance=0.025):
    """
    utils/observed_mz
    description:
            A generalized version of utils/monoisotopic_peak. Identifies the observed (i.e., feature) m/z peak and isotopologue type within the MS1 scan. This function iteratively compares the intensity of the target feature with peaks at m/z - 1 Da, m/z - 2 Da, etc. The peak with the highest intensity is assumed to be the monoisotopic peak.
    parameters:
            identified_peaks (list of tuples) -- extracted and centroided MS1 spectrum.
            target_mz (float) -- target precursor m/z value to compare with extracted peaks.
            tolerance (float) -- Da tolerance around the target m/z value to consider. Default is 0.025 Da.
    returns:
            (tuple) -- extracted (m/z, intensity) values of the identified monoisotopic peak.
    """

    # Filter peaks within the m/z tolerance range
    potential_peaks = [
        peak for peak in identified_peaks if abs(peak[0] - target_mz) <= mz_tolerance
    ]

    # If no peaks are found within tolerance, return the target m/z
    if not potential_peaks:

        return (target_mz, 0), 0

    # Assume the most intense peak is the monoisotopic peak
    monoisotopic_peak = max(potential_peaks, key=lambda x: x[1])

    # Initialize variables to hold the identified monoisotopic peak and its isotopologue type
    isotopologue_type = 0

    # Check for isotopologue peaks (M + 1, M + 2, ...) and compare intensities
    # Check up to M + 3
    for delta_mz in range(1, 4):
        isotopologue_mz = target_mz - delta_mz
        isotopologue_peak = next(
            (
                peak
                for peak in potential_peaks
                if abs(peak[0] - isotopologue_mz) <= mz_tolerance
            ),
            None,
        )

        # Update the monoisotopic peak if an isotopologue is more intense
        if isotopologue_peak and isotopologue_peak[1] > monoisotopic_peak[1]:
            monoisotopic_peak = isotopologue_peak
            isotopologue_peak = delta_mz

    return monoisotopic_peak, isotopologue_type
