## xulab_software/imms/expysomics

This Python package is intended to facilitate the integration of data independent (DIA) IM-MS/MS data into nontarget and suspect screening workflows. The current version of `expysomics` supports the detection and identification of quaternary ammonium compounds (QACs) and their phase I hepatic metabolites in human biological samples; the complete QAC reference database referred to herein is available at https://ccsbase.net/qac. Please note that the `dhrmasslynxapi` package containing the original Waters sdk files, copied into the `sdk\` directory, are required to utilize the `expysomics` package. The sdk files may be obtained from the lab NAS (under `lab_resources/masslynx_sdk_files/`):
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
* `process.py`
* `multigauss.py`
* `filter.py`
* `query.py`
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
8. `multigauss.process_data`
9. `filter.filter_data`

### Initialization
The `RawProcessor` object is initialized with the path to an Excel (.xlsx) spreadsheet, which contains the list of target m/z values and associated experimental metadata:

```python
from expysomics.process import RawProcessor

# initialize RawProcessor object 
data = RawProcessor("path/to/target/list.xlsx")
```
In the current version of `expysomics`, the target list **must** contain the following headers and information:

| file_name | mz    | sample_type | ccs_calibrant | gradient             | column_type |   
| ----------|:-----:|:-----------:|:-------------:|:--------------------:|:-----------:|
| file1.raw | 100.4 | control     | polyala       | develop_20min_5B_85B | Kinetex_C8  |
| file2.raw | 100.4 | sample      | polyala       | develop_20min_5B_85B | Kinetex_C8  |


1. **file_name**: path to the Waters .raw files to be processed
2. **mz**: list of all the target m/z values used to extract spectral features
3. **sample_type**: indicate whether the Waters .raw file being processed is an experimental "control" or a "sample" 
4. **ccs_calibrant**: indicate whether the type of CCS calibrant used is "polyala," "drugs," or "agilent"
5. **gradient**: indicate the type of LC gradient used
6. **column_type**: indicate the type of LC column used

A template spreadsheet is provided within this repository for convenience.

### Methods
`RawProcessor.extract`

This is the primary method for extracting data from .raw files for downstream analysis and database matching. 

#### Parameters
* `calibration_file` `(str)` - path to the CCS calibration file (described below)
* `ms1_function` `(float)` - function to use for MS1 scan
* `mobility_function` `(float)` - function containing mobility data
* `mz_tolerance` `(float)` - mass tolerance (in Da)

To convert extracted drift times into calibrated CCS values, the path to an Excel (.xlsx) spreadsheet with the following information must be indicated:
1. **observed calibrant m/z vaues**: list of the observed m/z values for a set of CCS calibrants
2. **reference CCS values**: literature CCS values of the calibrants used 

An example CCS calibrant spreadsheet is provided within this repository.

#### Returns
Excel (.xlsx) spreadsheet containing the following extracted data for each target m/z:

* monoisotopic_mz: m/z value of the monoisotopic peak (method adapted from MetaboAnnotatoR)
* dt: retention time-selected drift time (if available)
* ccs: calibrated CCS value (if drift time is available)
* rt: retention time 
* peak_area: peak area of the extracted LC ion chromatogram

If the input target list contains the paths to .raw files corresponding to a blank control, `RawProcessor.extract` returns an additional Excel (.xlsx) spreadsheet containing filtered spectral features (i.e., those that are unique to samples based on LC ion chromatogram peak area).

Finally, `RawProcessor.extract` creates two folders in the root directory (if they do not already exist) called "Extracted Chromatograms" and "Extracted Mobilograms" where .png images of the extracted, smoothed, and fitted data will be saved.

#### Example
```python
# extract data from .raw files indicated in target list
# function 1 contains MS1 data
# function 2 contains mobility data
# utilize a +/- 0.025 Da mass tolerance 
data = RawProcessor.extract("path/to/calibration/file.xlsx", 1, 2, 0.025)
```

## `expysomics.query`
This module defines the FeatureAnnotate object for matching extracted spectral features with potential database compounds. 

### Dependencies
1. `pandas`
2. `sqlite3`
3. `numpy`
4. `math`
5. `os`
6. `matplotlib`
7. `scipy`
8. `multigauss.process_data`
9. `dhrmasslynxapi.reader`

### Initialization
The `FeatureAnnotate` object is initialized with the paths to a spectral feature list (i.e., the output of `RawProcessor.extract`) and a database (.db) file containing reference experimental values. **If MS/MS spectral matching is also desired, the path to a spectral database (.db) file must also be indicated.**

```python
from expysomics.query import FeatureAnnotate

# initialize FeatureAnnotate object 
matches = FeatureAnnotate("path/to/feature/list.xlsx", "path/to/reference/db.db" "optional/path/to/spectral/db.db")
```

The current version of `expysomics` is set up to accept the output file of `RawProcessor.extract` without further modifications. To utilize this module with custom datasets, the spectral feature list **must** contain the following headers and information: 

1. **file_name**: path to the Waters .raw files to be processed
2. **mz**: list of all the target m/z values used to extract spectral features
3. **sample_type**: indicate whether the Waters .raw file being processed is an experimental "control" or a "sample" 
	* note that only spectral features identified as "sample" will be matched
4. **ccs_calibrant**: indicate whether the type of CCS calibrant used is "polyala," "drugs," or "agilent"
5. **gradient**: indicate the type of LC gradient used
6. **column_type**: indicate the type of LC column used
7. **dt**: drift time of the spectral feature
8. **ccs** calibrated CCS value of the spectral feature
9. **rt** retention time of the spectral feature
10. **peak_area** peak area of the extracted LC ion chromatogram

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
This is the primary method for matching extracted spectral features with potential compounds within the QAC reference database based on m/z, retention time, CCS, and MS/MS spectral matching selection criteria.

`FeatureAnnotate.match_msms` annotates extracted spectral features with potential matches based on m/z, retention time, and CCS selection criteria. The reference MS/MS spectrum for each potential match (if any) is first retrieved from the spectral database indicated during the initialization of the FeatureAnnotate object. Then, all ion fragmentation (AIF) spectra are **(1)** extracted from the .raw data file around the feature's retention and drift time windows, **(2)** normalized, **(3)** deconvoluted (by attempting to align fragment ion drift times with the precursor ion drift time), and then **(4)** compared to reference spectra using the cosine similarity score. `FeatureAnnotate.match_msms` also returns a composite score, adapted from MetaboAnnotatoR (Ebbels et al., 2022), calculated by summing the **(1)** reciprocal of the m/z error between the target and potential match parent, **(2)** CCS difference between the target and potential match parent, and **(3)** similarity of the fragmentation patterns. Similar to MetaboAnnotatoR, `expysomics` downweighs, by a factor of 2, fragment ions that are not drift time-aligned with their precursor ion but still present in the extracted AIF spectrum. 

#### Parameters
* `ms1_function` `(float)` - function to use for MS1 scan
* `mobility_function` `(float)` - function containing mobility data
* `database_type` `(str)` - indicate use of "experimental" or "theoretical" spectral database
* `mz_tolerance` `(float)` - mass tolerance for database query and AIF spectrum processing (in Da)
* `rt_tolerance` `(float)` - retention time tolerance for database query (in min)
* `ccs_tolerance` `(float)` - CCS tolerance for database query (in %)

#### Returns
* Excel (.xlsx) spreadsheet with potential matches (if any), ranked by their composite score (if multiple potential matches for a given specrtal feature)
* mirror plot (.png) displaying the reference fragmentation spectrum and deconvoluted AIF extracted spectrum

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
