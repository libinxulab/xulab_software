"""
    toxccs/utils/threshold.py

    Ryan Nguyen (ryan97@uw.edu)
    12/20/24

    description:
            Module with utility for evaluating Gaussian-fitted extracted ion mobilograms (EIMs).
"""


def fwhm_threshold(C, dt_i, fwhm_thresholds=(0.05, 2.5), intensity_threshold=500):
    """
    toxccs/utils/threshold.py/fwhm_threshold
    description:
            Checks if the full-width at half maximum (FWHM) and intensity of the fitted EIM meet threshold values. Default minimum and maximum FWHM thresholds are 0.05 and 2.5 ms, respectively. Default intensity threshold is 500.
    parameters:
            C (float) -- standard deviation of the Gaussian curve.
            dt_i (array) -- intensity values.
            fwhm_thresholds (tuple) -- minimum and maximum acceptable FWHM values.
            intensity_threshold (float) -- minimum acceptable intensity.
    returns:
            (bool) -- True if the FWHM and intensity meet the thresholds. False otherwise.
    """

    fwhm = C * 2.355
    max_intensity = max(dt_i)

    return (
        fwhm_thresholds[0] < fwhm < fwhm_thresholds[1]
        and max_intensity > intensity_threshold
    )
