## toxccs

This Python package is intended to automate the extraction and processing of LC-IM-MS/MS data for building multi-dimensional reference databases. The current version of `toxccs` interfaces directly with Waters .raw data and extracts both chromatographic and mobility information for manual review. For Please note that the `dhrmasslynxapi` package containing the original Waters sdk files, copied into the `sdk\` directory, are required to utilize the `toxccs` package and its sub-packages. The sdk files may be obtained from the lab NAS (under `lab_resources/masslynx_sdk_files/`):

* `cdt.dll`
* `MassLynxLockMassProcessor.py`
* `MassLynxRaw.dll`
* `MassLynxRaw.lib`
* `MassLynxRawChromatogramReader.py`
* `MassLynxRawDefs.py`
* `MassLynxRawInfoReader.py`
* `MassLynxRawReader.py`
* `MassLynxRawScanReader.py`
* `MassLynxSampleList.py`

The `dhrmasslynxapi` is located under xulab_software/imms and may be copied directly into the root directory. 

## `Directory Structure`
The root directory should contain the following Python files and folders:
* `dhrmasslynxapi` 
    * ensure that all files listed above are copied into the sdk\ directory
 * `toxccs`
* `Waters .raw files`
	* it is recommended to copy all required .raw files into the root directory from an external hard drive
* `Input .xlsx files` (described below)

## `toxccs.extract`
This module defines the RawProcessor object for extracting chromatographic and mobility information from Waters .raw files. Note that this module is intended to handle simple ion mobility profiles, where all peaks are fit to a single Gaussian function. If your mobility data are displaying more complex behavior, such as multiple peaks arising from protomers or constitutional isomers, it is highly recommended to use the `toxccs.multi_extract` module defined below.

### Dependencies
1. `pandas`
2. `numpy`
3. `scipy == 1.10.0`
5. `matplotlib`
6. `os`
7. `dhymasslynxapi.reader` and `dhrmasslynxapi.ccs_calibration`

### Initialization
The `RawProcessor` object is initialized with the path to an Excel file (.xlsx), which contains the list of target m/z values and associated experimental metadata:

```python
from toxccs.extract import RawProcessor

# initialize RawProcessor object 
data = RawProcessor("path/to/target/list.xlsx")
```
In the current version of `toxccs`, the target list should contain the following headers and information:

| File Name | Compound Name | SMILES                | Adduct | Target m/z | Sample Type | CCS Calibrant | Gradient     | Column Type |  
| ----------|:-------------:|:---------------------:|:------:|:----------:|:-----------:|:-------------:|:------------:|:-----------:|
| file1.raw | Acetaminophen | CC(=O)NC1=CC=C(C=C1)O | [M+H]+ | 152.0708   | sample      | polyala       | develop_100B | Kinetex_C18 |
| file2.raw | Acetaminophen | CC(=O)NC1=CC=C(C=C1)O | [M-H]- | 150.0563   | sample      | polyala       | develop_100B | Kinetex_C18 |

1. **File Name**: (REQUIRED) path to the Waters .raw files to be processed
2. **Compound Name**: (REQUIRED) name of the target compound to be processed 
3. **SMILES** (OPTIONAL) SMILES string of the target compound to be processed
4. **Adduct** (REQUIRED) adduct type (i.e., [M+H]+, [M-H]-, etc.) of the target compound to be processed
5. **Target m/z**: (REQUIRED) theoretical m/z value of the target compound to be processed
6. **Sample Type**: (REQUIRED) indicate the sample type being processed (i.e., experimental "control," "sample," or  internal standard "IS")
7. **CCS Calibrant**: (REQUIRED) indicate the type of CCS calibrant used (i.e., "polyala," "drugs," "agilent," etc.)
8. **Gradient**: (REQUIRED) indicate the type of LC gradient used
9. **Columm Type**: (REQUIRED) indicate the type of LC column used

A template spreadsheet for creating the target m/z list is provided within this repository for convenience.

In addition to the required target m/z list, the `RawProcessor` object may also be initialized with the path to an optional Excel (.xlsx) spreadsheet containing reference retention times and CCS values for internal standards or quality control compounds:

```python
from toxccs.extract import RawProcessor

# initialize RawProcessor object with additional reference retention time and CCS list for quality control purposes
data = RawProcessor("path/to/target/list.xlsx", "optional/path/to/reference/list.xlsx")
```
If included, the reference list **must** contain the following headers:

| Exact m/z | Internal Standard | Gradient     | Reference Retention Time (min) | Reference CCS (Å²) |        
| ----------|:-----------------:|:------------:|:------------------------------:|:------------------:|
| 107.1234  | IS_1              | develop_100B | 10.02                          | 179.4              | 
| 207.1234  | IS_2              | develop_100B | 11.57                          | 189.76             | 

### Methods
`RawProcessor.extract`
This is the primary method for extracting and processing LC-IM-MS/MS data.

#### Parameters
* `calibration_file` `(str)` - path to the CCS calibration file (described below)
* `ms1_function` `(float)` - function containing MS1-level data
* `mobility_function` `(float)` - function containing mobility data
* `mz_tolerance` `(float)` - mass tolerance (in Da)

To convert extracted drift times into calibrated CCS values, the path to an Excel (.xlsx) spreadsheet with the following information **must** be indicated:
1. **observed calibrant m/z vaues**: list of the observed m/z values for a set of CCS calibrants
2. **reference CCS values**: literature CCS values of the calibrants used 
3. **path to the Waters .raw file**: indicated by the **sheet name**

A template spreadsheet for creating the CCS calibrant list is provided within this repository for convenience. 

#### Returns
Excel (.xlsx) spreadsheet containing the following extracted data for each target m/z:

* Observed m/z: m/z value of the identified monoisotopic peak (method adapted from MetaboAnnotatoR; https://pubs.acs.org/doi/10.1021/acs.analchem.1c03032)
* Observed Drift Time (ms): retention time-selected drift time (if available)
* Observed CCS (Å²): calibrated CCS value (if drift time is available)
* Observed Retention Time (min): retention time 
* EIC Peak Intensity: maximum peak height of the extracted LC ion chromatogram
* EIC Peak Area: peak area of the extracted LC ion chromatogram

Finally, `RawProcessor.extract` creates a folder in the root directory (if it does not already exist) called "Extracted Data," where images of the extracted, smoothed, and Gaussian-fitted MS1, chromatographic, and mobility data are saved to. 

#### Example
```python
# Specify path to input target list
# Optional path to input internal standard list is not shown
# The lists MUST be Excel files (.xlsx)
# Refer to example target list files for template
data = RawProcessor("test.xlsx")

# Execute main extraction sequence for non-complex mobility data
# NOTE: All extracted ion mobilograms are handled as single Gaussian peaks 
# Speciy path to the input CCS calibration file
# Refere to example calibration file for template
# Specify MS1 function, mobility function, and m/z tolerance in Da in that order
data.extract("calibration_test.xlsx", ms1_function=0, mobility_function=0, mz_tolerance=0.025)
```

## `toxccs.multi_extract`
This module defines the MultiRawProcessor object for extracting complex chromatographic and mobility information from Waters .raw files. Note that this module is recommended for most cases, especially where the mobility data have not been pre-processed or manually reviewed beforehand. 

### Dependencies
1. `pandas`
2. `numpy`
3. `scipy == 1.10.0`
4. `matplotlib`
5. `os`
6. `dhymasslynxapi.reader` and `dhrmasslynxapi.ccs_calibration`

### Initialization
The `MultiRawProcessor` object is also initialized with the path to an Excel file (.xlsx) containing the list of target m/z values and associated experimental metadata. **This input file should have the same headers as the Exel spreadsheet described above for initializing the `RawProcessor` object:**

```python
from toxccs.multi_extract import MultiRawProcessor

# initialize MultiRawProcessor object 
data = MultiRawProcessor("path/to/target/list.xlsx")
```

### Methods
`MultiRawProcessor.multi_extract`

This is the primary method for extracting and processing complex LC-IM-MS/MS data.

#### Parameters
* `calibration_file` `(str)` - path to the CCS calibration file (described below)
* `ms1_function` `(float)` - function containing MS1-level data
* `mobility_function` `(float)` - function containing mobility data
* `mz_tolerance` `(float)` - mass tolerance (in Da)

To convert extracted drift times into calibrated CCS values, the path to an Excel (.xlsx) spreadsheet with the following information **must** be indicated as described above for `toxccs.extract`.

#### Returns
Excel (.xlsx) spreadsheet containing the following extracted data for each target m/z:

* Observed m/z: m/z value of the identified monoisotopic peak (method adapted from MetaboAnnotatoR; https://pubs.acs.org/doi/10.1021/acs.analchem.1c03032)
* Observed Drift Time (ms): retention time-selected drift time (if available)
* Observed CCS (Å²): calibrated CCS value (if drift time is available)
* Observed Retention Time (min): retention time 
* EIC Peak Intensity: maximum peak height of the extracted LC ion chromatogram
* EIC Peak Area: peak area of the extracted LC ion chromatogram
* LC FWHM Flag: flag indicating the extracted LC ion chromaogram exceeds the full-width at half maximum (FWHM) threshold; intended to be used for post-processing database review
* EIM FWHM (ms): full-width at half maximum of the extracted ion mobilogram
* Two-Peak EIM Resolution: calculated resolution between two neigboring extracted ion mobilograms; calculated value is always placed in the cell first peak of the two
 
Finally, `MultiRawProcessor.multi_extract` creates a folder in the root directory (if it does not already exist) called "Extracted Reprocessed Data," where images of the extracted, smoothed, and Gaussian-fitted MS1, chromatographic, and mobility data are saved to. 

#### Example
```python
# Specify path to input target list
# The list MUST be an Excel file (.xlsx)
# Refer to example target list file for template
data = MultiRawProcessor("extract_test.xlsx")

# Execute main extraction sequence for complex mobility data
# Specify MS1 function, mobility function, and m/z tolerance in Da in that order
data.multi_extract("calibration_test.xlsx", ms1_function=0, mobility_function=0, mz_tolerance=0.025)
```

## `toxccs.dda_extract`
This module defines the DDARawProcessor object for extracting MS/MS fragmentation spectra from Waters Fast-DDA .raw files. Note that the raw data **must** be centroided prior to implementation; for post-acqusition processing of continuum data, refer to the methods available in `toxccs.utils`.

### Dependencies
1. `pandas`
2. `numpy`
3. `scipy == 1.10.0`
4. `matplotlib`
5. `os`
6. `dhymasslynxapi.reader` 

### Initialization
The `DDARawProcessor` object is initialized with the path to an Excel file (.xlsx), which contains the list of target m/z values and associated experimental metadata:

```python
from toxccs.dda_extract import DDARawProcessor

# initialize DDARawProcessor object 
data = DDARawProcessor("path/to/target/list.xlsx")
```
In the current version of `toxccs`, the target list **must** contain the following headers and information:

1. **File Name**: (REQUIRED) path to the Waters .raw files to be processed
2. **Commpound Name** (REQUIRED) name of the target compound to be processed
3. **SMILES** (OPTIONAL) SMILES string of the target compound to be processed
4. **Adduct** (REQUIRED) adduct type (i.e., [M+H]+, [M-H]-, etc.) of the target compound to be processed
5. **Target m/z**: (REQUIRED) theoretical m/z value of the target compound to be processed
6. **Sample Type**: (REQUIRED) indicate the sample type being processed (i.e., experimental "control," "sample," or  internal standard "IS")
7. **Gradient**: (REQUIRED) indicate the type of LC gradient used
8. **Columm Type**: (REQUIRED) indicate the type of LC column used

A template spreadsheet for creating the target m/z list is provided within this repository for convenience.

### Methods
`RawProcessor.dda_extract`

This is the primary method for extracting MS/MS fragmentation spectra from Waters Fast-DDA .raw files.

#### Parameters
* `mz_tolerance` `(float)` - mass tolerance (in Da)
* `num_precursors` `(int)` - number of precursor ions selected for MS/MS fragmentation (in MS Method) 

Note that if the target list contains Fast-DDA files with different `num_precursors` specified in their MS Methods, **it is recommended that `num_precursors` be explicitly defined as the highest value in the entire dataset.** For example, if data are extracted from 2 unique .raw files with 3 and 15 precursors selected for fragmentation at a given MS/MS event, `num_precursors` should be set to 15. 

#### Returns
Individual Excel (.xlsx) spreadsheet containing the extracted MS/MS fragmentation spectrum **for each target compound:**

* m/z: list of m/z values
* Intensity: corresponding list of intensity values 

Excel (.xlsx) spreadsheet containing the following extracted data:

* Observed m/z: m/z value of the identified monoisotopic peak (method adapted from MetaboAnnotatoR; https://pubs.acs.org/doi/10.1021/acs.analchem.1c03032)
* Channel Number: value of the Channel from which the MS/MS fragmentation spectrum was extracted
* Mass error (ppm): mass error between the target m/z and observed m/z
* Flag: "FLAG" message for potentially false extracted MS/MS fragmentation spectra 

Finally, `DDARawProcessor.extract` creates a folder in the root directory (if it do not already exist) called "Extracted Data," where images of the chromatographic and MS/MS spectral data are saved to. 

#### Example
```python
# Specify path to input target list
# The list MUST be an Excel file (.xlsx)
# Refer to example target list file for template
data = DDARawProcessor("dda_test.xlsx")

# Execute main extraction sequence for MS/MS fragmentation spectra
# Data must be centroided and collected in Fast-DDA mode
# Specify m/z tolerance in Da and number of precursor ions selected for fragmentation
# Specify acquisition window (i.e., 0 to 1.8 min) for dynamic plotting window
data.dda_extract(mz_tolerance=0.01, num_precursors=5, x_min=0, x_max=1.8)
```





