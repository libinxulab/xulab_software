"""
test.py

Ryan Nguyen (ryan97@uw.edu)
12/24/24

description:
        Example script for utilizing sub-packages within toxccs. Simply uncomment the block(s) of code required for analysis and run the command: python3 test.py. 
"""

from toxccs.extract import RawProcessor
from toxccs.multi_extract import MultiRawProcessor
from toxccs.dda_extract import DDARawProcessor


# Specify path to input target list
# The list MUST be an Excel file (.xlsx)
# Refer to example target list file for template
# data = RawProcessor("extract_test.xlsx")

# Execute main sequence for non-complex mobility data
# NOTE: All extracted ion mobilograms are handled as single Gaussian peaks -- for handling complex data, use data.multi_extract below
# Speciy path to the input CCS calibration file
# Refere to example calibration file for template
# Specify MS1 function, mobility function, and m/z tolerance in Da
# data.extract("calibration_test.xlsx", ms1_function=0, mobility_function=0, mz_tolerance=0.025)

# Execute main sequence for complex mobility data
# data = MultiRawProcessor("test.xlsx")
# data.multi_extract("calibration_test.xlsx", ms1_function=0, mobility_function=0, mz_tolerance=0.025)

# Execute main sequence for MS/MS fragmentation spectra
# Data must be centroided and collected in Fast-DDA mode
# Specify m/z tolerance in Da and number of precursor ions selected for fragmentation
# Specify acquisition window (i.e., 0 to 1.8 min) for dynamic plotting window
# data = DDARawProcessor("dda_test.xlsx")
# data.dda_extract(mz_tolerance=0.025, num_precursors=5, x_min=0, x_max=1.8)
