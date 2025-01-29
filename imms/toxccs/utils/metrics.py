"""
toxccs/utils/metrics.py

Ryan Nguyen (ryan97@uw.edu)
12/21/24

description:
        Sub-package with utility for calculating mass spectrometric peak metrics (i.e, resolution and mass accuracy).
"""


def calculate_peak_resolution(dt_a, fwhm_a, dt_b, fwhm_b):
    """
    toxccs/utils/metrics.py/calculate_peak_resolution
    description:
        Calculates the resolution between adjacent IM peaks using the formula decribed by Dodds et al. (2017).
    parameters:
        dt_a (float or str) -- drift time of peak A.
        fwhm_a (float or str) -- FWHM value of peak A.
        dt_b (float or str) -- drift time of peak B.
        fwhm_b (float or str) -- FWHM value o peak B.
    returns:
        resolution_values (float) -- calculated resolution between the pair of peaks.
    """

    try:

        # Strip leading/trailing spaces from strings and handle empty strings
        dt_a = dt_a.strip() if isinstance(dt_a, str) else dt_a
        fwhm_a = fwhm_a.strip() if isinstance(fwhm_a, str) else fwhm_a
        dt_b = dt_b.strip() if isinstance(dt_b, str) else dt_b
        fwhm_b = fwhm_b.strip() if isinstance(fwhm_b, str) else fwhm_b

        # Ensure that all inputs are numbers and valid
        if any(v in [None, ""] for v in [dt_a, fwhm_a, dt_b, fwhm_b]):
            return None
        dt_a, fwhm_a, dt_b, fwhm_b = map(float, [dt_a, fwhm_a, dt_b, fwhm_b])

        # Convert values to floats
        dt_a = float(dt_a)
        fwhm_a = float(fwhm_a)
        dt_b = float(dt_b)
        fwhm_b = float(fwhm_b)

        # Avoid division by 0
        if (fwhm_a + fwhm_b) == 0:
            return None

        return 1.18 * abs(dt_b - dt_a) / (fwhm_a + fwhm_b)
    except (ValueError, TypeError) as e:
        print(f"Error calculating resolution: {e}")
        return None


def calculate_mass_error(theoretical_mz, observed_mz):
    """
    toxccs/utils/metrics.py/calculate_mass_error
    description:
            Calculates mass error in parts per million (ppm).
    parameters:
            theoretical_mz (float) -- theoretical m/z value.
            observed_mz (float) -- observed m/z value.
    returns:
            (float) -- mass error in ppm.
    """
    return ((observed_mz - theoretical_mz) / theoretical_mz) * 1e6
