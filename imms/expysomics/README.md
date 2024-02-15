## xulab_software/imms/expysomics

This Python package is intended to facilitate the efficient integration of multi-dimensional, data independent (DIA) IM-MS/MS data into suspect screening workflows. The current version of `expysomics` supports the detection and identification of quaternary ammonium compounds (QACs) and their phase I hepatic metabolites in human biological samples; the complete QAC reference database referred to herein is available at https://ccsbase.net/qac. Please note that the `dhrmasslynxapi` package containing the original Waters sdk files, copied into the `sdk\` directory, are required to utilize the `expysomics` package. The sdk files may be obtained from the lab NAS (under `lab_resources/masslynx_sdk_files/`):

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
 * `expysomics`
	 * `process.py`
	 * `multigauss.py`
	 * `multigauss.py`
* `qacs.db` 
* `qacs_experimental_msms.db`
* `qacs_theoretical_msms.db`
* `Waters .raw files`
* `Input .xlsx files` (described below)

## `expysomics.process`
This module defines the RawProcessor object for extracting data from .raw files.

### Dependencies
1. `pandas`
2. `numpy`
3. `scipy`
4. `math`
5. `matplotlib`
6. `os`
7. `dhymasslynxapi.reader` and `dhrmasslynxapi.ccs_calibration`
8. `multigauss.process_chromatogram`
9. `filter.filter_data`

### Initialization
The `RawProcessor` object is initialized with the path to an Excel (.xlsx) spreadsheet, which contains the list of target m/z values and associated experimental metadata:

```python
from expysomics.process import RawProcessor

# initialize RawProcessor object 
data = RawProcessor("path/to/target/list.xlsx")
```
In the current version of `expysomics`, the target list **must** contain the following headers and information:

| File Name | Target m/z | Sample Type | CCS Calibrant | Gradient             | Column Type |   
| ----------|:----------:|:-----------:|:-------------:|:--------------------:|:-----------:|
| file1.raw | 100.4      | control     | polyala       | develop_20min_5B_85B | Kinetex_C8  |
| file2.raw | 100.4      | sample      | polyala       | develop_20min_5B_85B | Kinetex_C8  |

1. **File Name**: path to the Waters .raw files to be processed
2. **Target m/z**: target m/z value used to extract spectral features
3. **Sample Type**: indicate the sample type being processed (i.e., experimental "control," "sample," or  internal standard "IS")
4. **CCS Calibrant**: indicate the type of CCS calibrant used (i.e., "polyala," "drugs," or "agilent")
5. **Gradient**: indicate the type of LC gradient used
6. **Columm Type**: indicate the type of LC column used

A template spreadsheet for creating the target m/z list is provided within this repository for convenience.

In addition to the required target m/z list, the `RawProcessor` object may also be initialized with the path to an optional Excel (.xlsx) spreadsheet containing reference retention times and CCS values for internal standards or quality control compounds:

```python
from expysomics.process import RawProcessor

# initialize RawProcessor object with additional reference retention time and CCS list for quality control purposes
data = RawProcessor("path/to/target/list.xlsx", "optional/path/to/reference/list.xlsx")
```
If included, the reference list **must** contain the following headers:

| Exact m/z | Internal Standard | Gradient              | Reference Retention Time (min) | Reference CCS (Å²)   |        
| ----------|:-----------------:|:---------------------:|:------------------------------:|:--------------------:|
| 283.3126  | d7-C10-BAC        | develop_20min_5B_85B  | 10.02                          | 179.4                | 
| 311.3439  | d7-C12-BAC        | develop_20min_5B_85B  | 11.57                          | 189.76               | 

A template spreadsheet for creating the reference list is provided within this repository for convenience. 

### Methods
`RawProcessor.extract`

This is the primary method for extracting data from .raw files for downstream analysis and database matching. 

#### Parameters
* `calibration_file` `(str)` - path to the CCS calibration file (described below)
* `ms1_function` `(float)` - function containing MS1-level data
* `mobility_function` `(float)` - function containing mobility data
* `mz_tolerance` `(float)` - mass tolerance (in Da)

To convert extracted drift times into calibrated CCS values, the path to an Excel (.xlsx) spreadsheet with the following information must be indicated:
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

If the `RawProcessor` object is initialized with the path to a reference list, the output spreadsheet will also flag internal standard (IS) features that deviate significantly (greater than 0.2 minutes and 3% for retention time and CCS, respectively) from the reference values. This functionality provides a simple method for assessing the quality and reproducibility of the chromatographic and mobility separations. The current version of `expysomics` does not support automatic retention time correction.

Additionally, if the input target list contains the paths to .raw files corresponding to a blank control, `RawProcessor.extract` returns an additional Excel (.xlsx) spreadsheet containing filtered spectral features (i.e., those that are unique to samples based on LC ion chromatogram peak area). An additional column will also be displayed that indicates the EIC peak area ratio between the sample and the corresponding control features, if detected. 

Finally, `RawProcessor.extract` creates two folders in the root directory (if they do not already exist) called "Extracted Chromatograms" and "Extracted Mobilograms" where .png images of the extracted, smoothed, and fitted data are saved to. 

#### Example
```python
# extract data from .raw files indicated in target list
# function 1 contains MS1 data
# function 2 contains mobility data
# utilize a +/- 0.025 Da mass tolerance to extract features
data.extract("path/to/calibration/file.xlsx", 1, 2, 0.025)
```

## `expysomics.query`
This module defines the FeatureAnnotate object for matching extracted spectral features with potential database compounds. 

### Dependencies
1. `pandas`
2. `sqlite3`
3. `numpy`
4. `math`
5. `os`
6. `re`
7. `matplotlib`
8. `scipy`
9. `multigauss.process_chromatogram`
10. `dhrmasslynxapi.reader`
11. `brainpy`

### Initialization
The `FeatureAnnotate` object is initialized with the paths to a spectral feature list (i.e., the output of `RawProcessor.extract`) and a database (.db) file containing reference experimental values. **If MS/MS spectral matching is also desired, the path to a spectral database (.db) file must also be indicated.**

```python
from expysomics.query import FeatureAnnotate

# initialize FeatureAnnotate object 
matches = FeatureAnnotate("path/to/feature/list.xlsx", "path/to/reference/db.db" "optional/path/to/spectral/db.db")
```

Note that the QAC reference and spectral databases are included within this repository.

The current version of `expysomics` is designed to accept the output file of `RawProcessor.extract` without further modifications. To utilize this module with custom datasets, the input spectral feature list **must** contain the following headers and information: 

1. **File Name**: path to the Waters .raw files to be processed
2. **Sample Type** indicate the sample type being processed (i.e., experimental "control," "sample," or  internal standard "IS")
3. **Gradient** indicate the type of LC gradient used
4. **Column Type** indicate the type of LC column used
5. **Target m/z** target m/z value used to extract spectral features
	* note that only spectral features identified as "sample" will be matched against the database
6. **Observed m/z** m/z value of the identified monoisotopic peak
7. **CCS Calibrant**: indicate the type of CCS calibrant used (i.e., "polyala," "drugs," or "agilent")
8. **Observed Drift Time (ms)**: drift time of the spectral feature
9. **Observed CCS (Å²)** calibrated CCS value of the spectral feature
10. **Observed Retention Time (min)** retention time of the spectral feature
11. **EIC Peak Intensity** maximum peak height of the extracted LC ion chromatogram
12. **EIC Peak Area** peak area of the extracted LC ion chromatogram
13. **EIC Peak Area Ratio** ratio, if available, between the peak area of the sample feature and the corresponding control feature (i.e., the factor by which the sample feature peak area is above the control)

### Methods
`FeatureAnnotate.match_features`
This is the primary method for matching extracted spectral features with potential database compounds based on various selection criteria (i.e., m/z, retention time, and CCS).

#### Parameters
* `mz_tolerance` `(float)` - mass tolerance for database query (in Da)
* `rt_tolerance` (`None` or `float`) - if not `None`, retention time tolerance for database query (in min) 
* `ccs_tolerance` (`None` or `float`) - if not `None`, CCS tolerance for database query (in %)

The feature list will be annotated at the highest criteria level indicated; i.e., if both `rt_tolerance` and `ccs_tolerance` are passed to `FeatureAnnotate.match_features`, potential matches will be fetched based on m/z, retention time, and CCS criteria using the indicated tolerances for each parameter. Note that `mz_tolerance` is always required. 

#### Returns
* Excel (.xlsx) spreadsheet containing the potential matches, if any, for the feature list.

#### Example
```python
# annotate feature list based on m/z, retention time, and CCS selection criteria
# utilize +/- 0.025 Da m/z tolerance
# and +/- 0.2 min retention time tolerance
# and +/- 3% CCS tolerance
matches.match_features(0.025, 0.2, 0.03)
```

`FeatureAnnotate.match_msms`
This is the primary method for matching extracted spectral features with potential compounds within the QAC reference database based on m/z, retention time, CCS, and MS/MS spectral matching selection criteria. **Note that the path to a spectral database (.db) file must be indicated during the initialization of `FeatureAnnotate` in order to utilize this function.**

`FeatureAnnotate.match_msms` annotates extracted spectral features with potential matches based on m/z, retention time, and CCS selection criteria. The reference MS/MS spectrum for each potential match (if any) is first retrieved from the spectral database indicated during the initialization of the FeatureAnnotate object. All ion fragmentation (AIF) spectra are then **(1)** extracted from the .raw data file around the feature's retention and drift time windows, **(2)** centroided, **(3)** deconvoluted (by attempting to align fragment ion drift times with the precursor ion drift time), and then **(4)** compared to reference spectra using the reverse dot product-based cosine similarity score. `FeatureAnnotate.match_msms` also returns a composite score, adapted from LPPTiger (https://www.nature.com/articles/s41598-017-15363-z), calculated by summing the **(1)** reciprocal of the m/z error between the target and potential match parent, **(2)** CCS difference between the target and potential match parent, **(3)** similarity of the fragmentation patterns, and **(4)** differences between the observed and theoretical isotopic distribution patterns in the MS1 scan. Similar to MetaboAnnotatoR, `expysomics` downweighs ( by a factor of 2) fragment ions that are not drift time-aligned with their precursor ion but still present in the extracted AIF spectrum. 

#### Parameters
* `ms1_function` `(float)` - function to use for MS1 scan
* `mobility_function` `(float)` - function containing mobility data
* `database_type` `(str)` - indicate use of "experimental" or "theoretical" spectral database
* `mz_tolerance` `(float)` - mass tolerance for database query and AIF spectrum processing (in Da)
* `rt_tolerance` `(float)` - retention time tolerance for database query (in min)
* `ccs_tolerance` `(float)` - CCS tolerance for database query (in %)

#### Returns
* Excel (.xlsx) spreadsheet with potential matches (if any), ranked by their composite score (if multiple potential matches for a given specrtal feature), and associated scores
* mirror plot (.png) displaying the reference fragmentation spectrum and deconvoluted AIF extracted spectrum
* image (.png) displaying the theoretical and observed isotopic distributions in the MS1 scan

#### Example
```python
# annotate feature list based on m/z, retention time, and CCS selection criteria
# function 1 contains MS1 data
# function 2 contains mobility data
# utilize experimental spectral database
# utilize +/- 0.025 Da m/z tolerance
# and +/- 0.2 min retention time tolerance
# and +/- 3% CCS tolerance
matches.match_msms(1, 2, "experimental", 0.025, 0.2, 0.03)
```
