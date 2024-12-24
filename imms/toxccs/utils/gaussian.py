"""
toxccs/utils/gaussian.py

Ryan Nguyen (ryan97@uw.edu)
12/20/24

description:
        Module with utilities for fitting raw mass spectrometry data to Gaussian functions. 
"""

from scipy.signal import convolve, gaussian, find_peaks
from scipy.optimize import curve_fit
import numpy as np


def multi_gaussian_fixed_mu(x, mu_values, *params):
    """
    toxccs/utils/gaussian.py/multi_gaussian_fixed_mu
    description:
            Computes the sum of multiple Gaussian functions with fixed mean (mu) values for each peak. This function is intended to primarily be used for fitting chromatographic and complex mobility data, where the mean positions of the peaks are known but their amplitudes and widths (i.e., standard deviations) are to be determined.
    parameters:
            x (numpy.ndarray) -- independent variable (e.g., rt or m/z) where the function is being evaluated.
            mu_values (numpy.ndarray) -- array of fixed mean (mu) values for the Gaussian functions. Each value corresponds to the mean position of a peak in the data.
            *params (list) -- parameters for the Gaussian functions, with each peak's amplitude and standard deviation.
    returns:
            y (numpy.ndarray) -- the evaluated sum of Gaussian functions at each point in x.
    """

    y = np.zeros_like(x)
    for i in range(0, len(params), 2):

        # Calculate the sum of multiple Gaussian functions
        A = params[i]
        sigma = params[i + 1]
        mu = mu_values[int(i / 2)] if len(mu_values) > int(i / 2) else 0
        y += A * np.exp(-((x - mu) ** 2) / (2 * sigma**2))

    return y


def gaussian_fit(x, A, B, C):
    """
    toxccs/utils/gaussian.py/gaussian_fit
    description:
            Single Gaussian function used for fitting extracted peaks.
    parameters:
            x (array) -- independent variable (e.g., t or m/z) where the function is being evaluated.
            A (float) -- amplitude of the Gaussian curve.
            B (float) -- mean value of the Gaussian curve.
            C (float) -- standard deviation of the Gaussian curve.
    returns:
            (array) -- evaluated Gaussian function at each point in x.
    """

    if abs(C) < 0.01:
        C = 0.01

    return A * np.exp(-((x - B) ** 2) / (2 * C**2))


def gaussian_smooth_pick(
    mz_array, intensity_array, window_len=1, std=0.1, prominence=0.001
):
    """
    toxccs/utils/gaussian.py/gaussian_smooth_pick
    description:
            Convolves the intensity array with a Gaussian window to smooth the data and identifies local maxima (i.e., peaks) based on a defined prominence. This function is analogous to centroiding the profile data, as outlined in MSnbase.
    parameters:
            mz_array (numpy.ndarray) -- array of m/z values.
            intensity_array (numpy.ndarray) -- array of intensity values.
            window_len (int) -- length of the Gaussian window used for smoothing. Default is 1.
            std (float) -- standard deviation of the Gaussian window, controlling the degree of smoothing. Default is 0.1.
            prominence (float) -- minimum prominence of peaks to be identified. Default is 0.001.
    returns:
            (tuple) -- a tuple containing identified peaks.
    """

    window = gaussian(window_len, std=std)
    smoothed_intensity = convolve(intensity_array, window / window.sum(), mode="same")

    # Identify peaks within the smoothed data
    peaks, _ = find_peaks(smoothed_intensity, prominence=prominence)
    identified_peaks = [(mz_array[p], smoothed_intensity[p]) for p in peaks]

    return identified_peaks, smoothed_intensity


def process_chromatogram(
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
    toxccs/utils/gaussian.py/process_chromatogram
    description:
        Processes raw data by smoothing time- or m/z-intensity pairs using Gaussian convolution, identifying peaks, and fitting the smoothed data to a multi-Gaussian function. Optionally generates and saves output .png images. Specifically written to handle input from single, non-complex mobilogram profiles,
    parameters:
        rt (list or numpy.ndarray) -- array of retention times.
        rt_i (list or numpy.ndarray) -- corresponding intensity values for each retention time.
        file_name (str) -- name of the .raw file from which the data is extracted.
        mz (float) -- target m/z value used to extract the raw data.
        sample_type (str) -- sample type (i.e., "control," "sample," or "IS").
        generation_images (bool, optional) -- flag to determine whether to generate and save chromatogram images. Default is false.
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
            while new_start_idx > 0 and rt_i[new_start_idx - 1] <= rt_i[new_start_idx]:
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
                lambda x, *p: multi_gaussian_fixed_mu(x, mu_values, *p),
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


def process_multi_chromatogram(
    rt, rt_i, prominence=200, distance=1, apply_fifty_percent_rule=True
):
    """
    toxccs/utils/gaussian.py/process_multi_chromatogram
    description:
        Processes raw data by smoothing time or m/z-intensity pairs using Gaussian convolution, identifying peaks, and
        fitting the smoothed data to a multi-Gaussian function. Specifically written to handle input from complex mobilogram profiles.
    parameters:
        rt (list or numpy.ndarray) -- array of retention times.
        rt_i (list or numpy.ndarray) -- corresponding intensity values for each retention time.
        prominence (float) -- threshold for how much the peak stands out from the surrounding baseline of the signal, defined as the vertical distance between the peak and its contour line. Default is 200.
        distance (float) -- threshold for minimum horizontal distance between neighboring peaks. Default (and minimum) is 1.
    returns:
        (tuple of (list of int, list of float, list of float, list of tuples))
            -- peak_indices -- indices of the identified peaks in the smoothed data.
            -- rt_list -- list of peaks corresponding to the identified peaks, rounded to 2 decimal places.
            -- areas -- calculated areas under each peak based on the multi-Gaussian fitting.
            -- peak_ranges -- start and end indices for each identified peak.
    """

    popt_multi_gaussian = []
    valid_popt_multi_gaussian = []
    window_len = 51
    window = gaussian(window_len, std=1)  # Adjust to change degree of smoothing
    smoothed_intensity = convolve(rt_i, window / window.sum(), mode="same")

    """# Initialize default values
        valid_peak_indices = []
        rt_list = []
        areas = []
        peak_ranges = []
        mu_values = []
        fwhm_flags = []
        fwhm_values = []"""

    # Identify peaks in the smoothed data
    # Default thresholds for prominence and distance are 100 and 1, respectively
    # These threshold values can be adjusted, and shoudld be kept relatively lenient, since identified peaks are subject to subsequent filtering based on the 10% valley rule
    peak_indices, _ = find_peaks(
        smoothed_intensity,
        prominence=prominence,
        distance=distance,
    )

    """if peak_indices.size == 0:
            return (
                valid_peak_indices,
                rt_list,
                areas,
                peak_ranges,
                valid_popt_multi_gaussian,
                mu_values,
                fwhm_flags,
                fwhm_values,
            )"""

    # Identify peak indices
    peak_ranges = []
    fwhm_flags = []
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
            while new_start_idx > 0 and rt_i[new_start_idx - 1] <= rt_i[new_start_idx]:
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
    rt_list = [str(round(float(rt[j]), 2)) for j in peak_indices]

    # Extract mu (mean) values for the fixed Gaussian functions from the smoothed data
    mu_values = rt[peak_indices] if peak_indices.size > 0 else np.array([])

    # Initial guesses for amplitude (A) and sigma from the smoothed data
    A_guesses = (
        smoothed_intensity[peak_indices] if peak_indices.size > 0 else np.array([])
    )
    sigma_guesses = [0.05] * len(peak_indices)
    p0 = [val for sublist in zip(A_guesses, sigma_guesses) for val in sublist]

    # Attempt to fit the smoothed data to a multi-Gaussian function
    try:
        if peak_indices.size > 0:
            popt_multi_gaussian, _ = curve_fit(
                lambda x, *p: multi_gaussian_fixed_mu(x, mu_values, *p),
                rt,
                rt_i,
                p0=p0,
                maxfev=20000,
            )

            # Apply 10% valley rule for identified peaks, if specified
            valid_peak_indices = peak_indices
            if apply_fifty_percent_rule:
                valid_peak_indices, valid_param_indices = fifty_percent_rule(
                    mu_values, popt_multi_gaussian, smoothed_intensity
                )
            else:
                valid_param_indices = list(range(len(popt_multi_gaussian)))
                valid_peak_indices = list(range(len(popt_multi_gaussian) // 2))

            # Map valid indices back to the original peak indices
            valid_peak_indices = [
                peak_indices[i] for i in valid_peak_indices if i < len(peak_indices)
            ]

            # Update mu_values based on valid peaks
            valid_mu_values = rt[valid_peak_indices]

            # Ensure correct mapping for popt_multi_gaussian
            valid_popt_multi_gaussian = [
                popt_multi_gaussian[i] for i in valid_param_indices
            ]

            # Calculate areas under the Gaussian curves
            areas = []
            fwhm_values = []
            for i in range(0, len(valid_popt_multi_gaussian), 2):
                A = valid_popt_multi_gaussian[i]
                sigma = valid_popt_multi_gaussian[i + 1]
                fwhm = sigma * 2.355
                mu = valid_mu_values[int(i / 2)]
                area = A * sigma * np.sqrt(2 * np.pi)  # Gaussian integral
                areas.append(area)
                fwhm_values.append(fwhm)

                # Handle cases where extracted LC ion chromatograms are wide, asymmetric, or otherwise exhibit significant deviation from ideal Gaussian behavior
                # This could be the case for analytes that elute over the entire chromatographic range due to mismatched LC conditions
                # When the FWHM of the peak is greater than 0.1 min, use the 5-95% bounded range of the fitted Gaussian curve for subsequent mobilogram extraction
                # This ensures that extracted mobilograms are not incorrectly truncated due to poor LC peak shape
                if fwhm > 0.1:
                    lower_bound = mu - 1.645 * sigma
                    upper_bound = mu + 1.645 * sigma
                    start_idx = np.argmin(np.abs(rt - lower_bound))
                    end_idx = np.argmin(np.abs(rt - upper_bound))

                    # Add FWHM flag for downstream identification
                    fwhm_flags.append("FWHM > 0.1 min")
                else:

                    # Revert to original retention time bounds when the extracted LC ion chromatograms are narrow and sharp
                    start_idx, end_idx = peak_ranges[int(i / 2)]
                    fwhm_flags.append("")
                peak_ranges[int(i / 2)] = (start_idx, end_idx)
            mu_values = valid_mu_values
        else:
            (
                valid_peak_indices,
                rt_list,
                areas,
                peak_ranges,
                valid_popt_multi_gaussian,
                mu_values,
                fwhm_flags,
                fwhm_values,
            ) = ([], [], [], [], [], [], [], [])
    except Exception as e:
        print(f"Exception occurred: {e}")
        areas = [None] * len(peak_indices)
        (
            valid_peak_indices,
            rt_list,
            peak_ranges,
            valid_popt_multi_gaussian,
            mu_values,
            fwhm_flags,
            fwhm_values,
        ) = ([], [], [], [], [], [], [])
        popt_multi_gaussian = []
    return (
        valid_peak_indices,
        rt_list,
        areas,
        peak_ranges,
        valid_popt_multi_gaussian,
        mu_values,
        fwhm_flags,
        fwhm_values,
    )


def peak_fit(t, dt_i, p0="guess"):
    """
    toxccs/utils/gaussian.py/peak_fit
    description:
            Fits a single peak in the extracted ion mobilogram (EIM) profile to a Gaussian curve.
    parameters:
            t (list or numpy.ndarray) -- drift times.
            dt_i (list or numpy.ndarray) -- intensity values.
            p0 (tuple or str) -- initial guesses for the Gaussian parameters.
    returns:
            (array) -- optimal values for the Gaussian parameters.
    """

    if p0 == "guess":
        p0 = (max(dt_i), t[np.argmax(dt_i)], 0.5)
    opt, cov = curve_fit(gaussian_fit, t, dt_i, maxfev=5000, p0=p0)

    return opt


def fifty_percent_rule(mu_values, popt_multi_gaussian, smoothed_intensity):
    """
    toxccs/utils/gaussian.py/apply_fifty_percent_rule
    description:
        Applies the 50% valley rule to idenify valid peaks in multi-Gaussian fitted data. The rule ensures that a peak is considered valid of the valley between it and adjacent peaks is sufficiently deep. Specifcally, the valid should be lower than the the higher peak's 10% height and heigher than the lower peak's 10% height. Peaks with a valley lower than both adjacent peaks' 10% height are also considered valid to ensure that isolated Gaussian peaks (i.e, those that are completely resolved from their neighboring peaks) are also considered valid. This filter ensures the accurate identification of peaks in cases of closely spaced or partially overlapping peaks.
    parameters:
        mu_values (numpy.ndarray) -- array of fixed mean (mu) values for the Gaussian functions.
        popt_multi_gaussian (list) -- optimized parameters for the multi-Gaussian fit, containing amplitudes and standard deviations of the Gaussian peaks.
        smoothed_intensity (numpy.ndarray) -- intensity values of the smoothed data.
    returns:
        valid_peak_indices (list) -- indices of the peaks that passed the 10% valley rule.
        valid_param_indices (list) -- indices of the corresponding Gaussian parameters for the valid peaks.
    """

    valid_peak_indices = []
    valid_param_indices = []
    num_peaks = len(mu_values)

    # Retrieve the full height for each Gaussian peak
    full_heights = [
        popt_multi_gaussian[i] for i in range(0, len(popt_multi_gaussian), 2)
    ]

    # Calculate the 10% height for each Gaussian peak
    fifty_percent_heights = [
        popt_multi_gaussian[i] * 0.5 for i in range(0, len(popt_multi_gaussian), 2)
    ]

    # Start from 1 to compare with the previous peak
    for i in range(1, num_peaks):
        prev_peak_idx = i - 1
        curr_peak_idx = i
        prev_peak_height = fifty_percent_heights[prev_peak_idx]
        curr_peak_height = fifty_percent_heights[curr_peak_idx]

        # Identify the valley between the two peaks
        # The valley is defined as the minimum intensity value in the smoothed data between the two peaks in question
        valley_min_idx = np.argmin(
            smoothed_intensity[
                int(mu_values[prev_peak_idx]) : int(mu_values[curr_peak_idx]) + 1
            ]
        )
        valley_min_idx += int(mu_values[prev_peak_idx])
        valley_height = smoothed_intensity[valley_min_idx]

        # Check if the valley is valid
        # The valley height must be lower than the higher of the two peaks' 10% heights
        # The valley height must also be higher than the lower of the two peaks' 10% heights
        # If the valley height is lower than both peaks' 10% heights, it is also considered valid (as in the case of fully resolved peaks)
        # Definitions adapted from Urban et al. (2014)
        if (
            (valley_height < max(prev_peak_height, curr_peak_height))
            and (valley_height > min(prev_peak_height, curr_peak_height))
            or (valley_height < prev_peak_height and valley_height < curr_peak_height)
        ):
            valid_peak_indices.append(prev_peak_idx)
            valid_peak_indices.append(curr_peak_idx)

            # Add both peaks as valid
            valid_param_indices.extend([2 * prev_peak_idx, 2 * prev_peak_idx + 1])
            valid_param_indices.extend([2 * curr_peak_idx, 2 * curr_peak_idx + 1])

    # Check for isolated peaks (those without neighboring peaks within a defined range)
    # May need to tweak this logic
    isolated_peak_indices = [
        i
        for i in range(num_peaks)
        if i == 0
        or i == num_peaks - 1
        or (i - 1 not in valid_peak_indices and i + 1 not in valid_peak_indices)
    ]

    # Filter isolated peaks by height threshold
    for idx in isolated_peak_indices:
        if full_heights[idx] >= 1000:
            valid_peak_indices.append(idx)
            valid_param_indices.extend([2 * idx, 2 * idx + 1])

    # Ensure unique peak indices
    valid_peak_indices = sorted(list(set(valid_peak_indices)))
    valid_param_indices = sorted(set(valid_param_indices))

    return valid_peak_indices, valid_param_indices
