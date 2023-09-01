## xulab_software/imms/qac_biomonitoring

These Python scripts are intended to facilitate the integration of IM-MS/MS into QAC biomonitoring experimental workflows. The `dhrmasslynxapi` package containing the original Waters sdk files, copied into the `sdk\` directory, are required. The sdk files may be obtained from the lab NAS (under `lab_resources/masslynx_sdk_files/`):
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

## `Directory Structure`
The root directory should contain the following Python files and folders:
* `dhrmasslynxapi` 
* `process.py`
* `process_dt_filter.py`
* `multigauss.py`
* `filter.py`
* `query.py`
* `qacs.db` 
* `qacs_experimental_msms.db`
* `qacs_theoretical_msms.db`
* `Waters .raw files`
* `Input .xlsx files` (described below)

## `process.py`
This script extracts the following spectral features from batches of Waters .raw files based on a list of target m/z values: 
* drift times (dt)
* calibrated collision cross section (CCS) values
* retention times (rt)
* LC ion chromatogram peak areas

**Note that `process.py` is intended to be used for independently extracting rt and dt, solely based on a target m/z value; if the user wishes to extract dt values using an appropriate rt filter, then `process_rt_filter.py` should be utilized instead.**  

### Dependencies
1. `pandas`
2. `numpy`
3. `scipy`
4. `math`
5. `matplotlib`
6. `os`
7. `dhymasslynxapi.reader` and `dhrmasslynxapi.ccs_calibration`
8. `multigauss`
9. `filter`

### Usage
To run `process.py`, create an Excel .xlsx spreadsheet containing the following headers and information:
1. **file_name**: List of all the Waters .raw files to be processed
2. **mz**: List of all the target m/z values used to extract spectral features
3. **sample_type**: Indicate whether the Waters .raw file being processed is a "control" or a "sample" 
4. **ccs_calibrant**: Indicate whether the type of CCS calibrant used is "polyala," "drugs," or "agilent"
5. **gradient**: Indicate the type of LC gradient used
6. **column_type**: Indicate the type of LC column used

**Note that column headers must exactly match those indicated here.** A template spreadsheet is provided within this repository for convenience.

**Additionally, an Excel .xlsx spreadsheet containing the following CCS calibrant information is required to convert extracted dt values into calibrated CCS values:**
1. **observed calibrant m/z values**: List of the observed m/z values for a set of CCS calibrants
2. **Reference CCS values**: Literature CCS values of the calibrants used and indicated 

**Note that sheet names must correspond to the exact, full name of the Waters .raw file containing CCS calibrant data. It is therefore highly recommended that both required spreadsheets be located in the same directory as all Waters .raw files and Python scripts described herein.** An example CCS calibrant spreadsheet is provided within this repository for convenience.

Before executing the script, paths to both spreadsheets must be explicitly indicated within `process.py`; there is currently no user interface. For example, if the script is run from the same directory containing the two input spreadsheets:
```python
if __name == "__main__":
	# Enter path to the processing spreadsheet template containing file names
	input_file = "input_file.xlsx"

	# Enter path to the CCS calibrant spreadsheet
	calibraiton_file = "calibrant_file.xlsx"
```

### Returns
* Excel .xlsx spreadsheet containing all extracted spectral features
* Excel .xlsx spreadsheet containing filtered spectral features unique to sample files (relative to control files, see `filter.py` below)
* Folder ("Extracted Mobilograms") containing .png images of extracted and fitted mobilograms
* Folder ("Extracted Chromatograms") containing .png images of extracted (raw), smoothed, and fitted LC ion chromatograms
* CCS calibration curve (.png) displaying randomly distributed fit residuals

### Notes
The current version of `process.py` utilizes multiple hardset parameters for interfacing with Waters .raw files. These parameters include MS function number(*e.g.* LC, IM) and desired mass tolerance around the target m/z. In order to override these parameters, the following line of code within `process.py` can be directly manipulated:
```python
# Extract mobilogram from MS function 2 for target m/z +/- 0.01 Da
rdr = MassLynxReader(file_name)
mz = float(mz)
t, dt_i = rdr.get_chrom(2, float(mz), 0.01)

# Extract m/z-selected EIC from MS function 1 for target m/z +/- 0.01 Da
rt, rt_i = rdr.get_chrom(1, float(mz), 0.01)
```
Note that in the first case, where a mobilogram is being extracted from the Waters .raw file, **the indicated MS function must contain mobility data.** 

Additionally, `process.py` utilizes two hardset filters to remove low intensity and/or poorly shaped extracted mobilograms/ATDs (based on full width at half maximum, FWHM). These two parameters may be overriden by directly manipulating the following lines of code:
```python
# Apply shape (0.01 ms < FWHM < 3 ms) and intensity (> 1000) thresholds on extracted mobilogram
...
if not (fwhm < 3 and fwhm > 0.01 and max(dt_i) > 1000): 
	...
```
## `process_rt_filter.py`
This script returns the same spectral features as `process.py` but implements an appropriate rt filter for extacting dt values. The same dependencies and Excel spreadsheets as described above are required to utilize `process_rt_filter.py`. 

This script is primarily intended to be used if the sample being processed is known or expected to contain complex mixtures of target compounds (*i.e.* human samples) with similar or identical m/z values. For example, if a sample contains two (or more) compounds that have m/z values within the input m/z tolerance window, the extracted LC ion chromatogram will contain multiple spectral features. Each of these spectral features presumably has their own dt and corresponding CCS value. The implementation of an appropriate rt filter is therefore necessary to distinguish spectral features that contain similar or identical m/z values; `process_rt_filter.py` does so by utilizing the `get_filtered_chrom` function within `dhrmasslynxapi`:
```python
# Extract mobilogram from MS function 1 for target m/z +/- 0.025 Da AND rt +/- 0.1 min of identified peak in extracted LC ion chromatogram
...
	t, dt_i = rdr.get_filtered_chrom(1, float(mz), 0.025, rt_min=float(rt_value)-0.1, rt_max=float(rt_value)+0.1)
		...
```

## `filter.py`
The main function of `filter.py` is to remove spectral features present in both sample files and their corresponding controls. This script is imported as a module within `process.py` and returns a list of unique spectral features.

### Dependencies
1. `pandas`

### Returns
* Excel .xlsx spreadsheet containing filtered spectral features unique to sample files

### Notes
The current version of `filter.py` removes non-unique spectral features on the basis of peak area. When two extracted features are aligned by both m/z and rt, their peak areas are compared: if the sample peak area is less than 1.5 times greater than that of the corresponding control peak area, it is rejected and not added to the filtered output spreadsheet. This peak area threshold can be overriden by directly manipulating the following line of code within `filter.py`:
```python
# Check if sample peak area is 2 times greater than control peak area
...
if sample_peak_area > 2 * control_peak_area:
	...
```

## `multigauss.py`
The `multigauss.py` script provides functionality for processing and analyzing LC chromatography data. This script **(1)** smooths raw extracted ion chromatograms (EICs) using Gaussian convolution, **(2)** identifies peaks within the smoothed data, and 
**(3)** simultaneously fits the raw EIC to multiple Gaussian functions for accurate peak identification and peak area calculations.

### Dependencies
* `numpy`
* `pandas`
* `os`
* `matplotlib`
* `scipy`

### Returns
* Folder ("Extracted Chromatograms") containing .png images of extracted (raw), smoothed, and fitted LC chromatograms

### Notes
While `multigauss.py` is intended to be imported and utilized within `process.py` or `process_rt_filter.py`, it may also be directly executed by manually inputing raw chromatographic data:
```python
from multigauss import process_data

# Manually input raw rt and intensity data
rt = [...]
rt_i = [...]

# Specify file name for output image
# Specify target m/z used to extract LC chromatogram
# Specify sample type (i.e. control or sample)
file_name = "sample_file_name"
mz = 123.45
sample_type = "sample"

# Process the raw chromatographic data 
peak_indices, rt_list, areas = process_data(rt, rt_i, file_name, mz, sample_type)
```

This code will smooth and fit the raw chromatographic data and return **(1)** the indices of identified peaks, **(2)** their rt, and **(3)** their peak areas.

There are several parameters within the current  version of `multigauss.py` that are hardset when imported and utilized within `process.py` and 'process_rt_filter.py`. For instance, the `find_peaks` module utilized to identify peaks within the smoothed data may be further optimized based on the characteristics of the input chromatographic data:
```python
# Adjust parameters to identify peaks within the smoothed chromatographic data
peak_indices, _ = find_peaks(smoothed_intensity, prominence = "", distance = "", width = "", height = "")
```
Even more parameters may be specified and adjusted within find_peaks; the current version of `multigauss.py` only utilizes the `prominence` and `distance` parameters.

## `query.py`
This script interfaces with the QACs Database for the purpose of querying and annotating LC-IM-MS/MS spectral features. The current version of `query.py` retrieves potential compound matches from the QAC Database based on user-defined criteria levels, including **(1)** m/z, **(2)** m/z and rt, **(3)** m/z and ccs, **(4)** m/z, rt, and ccs, and **(5)** m/z, rt, ccs, and MS/MS cosine similarity to stored reference spectra. 

### Dependencies
* `pandas`
* `sqlite3`
* `numpy`
* `math`
* `os`
* `matplotlib`
* `dhrmasslynxapi.reader`
* `qacs.db`, `qacs_experimental_msms.db`, and `qacs_theoretical_msms.db` (within the root directory, see below)

### Usages
The current version of `query.py` is specifically designed to accept filtered lists of spectral features, *i.e.* the output spreadsheet produced by `process.py` and/or `process_rt_filter.py`. If a list is obtained by other methods, the input spreadsheet must be reformmated to have the following structure:
1. **file_name**: List of all the Waters .raw files to be processed
2. **mz**: List of all the target m/z values used to extract spectral features
3. **sample_type**: Indicate whether the Waters .raw file is a "control" or a "sample" 
4. **ccs_calibrant**: Indicate whether the type of CCS calibrant used is "polyala," "drugs," or "agilent"
5. **gradient**: Indicate the type of LC gradient used
6. **column_type**: Indicate the type of LC column used
7. **dt**: The m/z-selected dt
8. **ccs**: The calibrated CCS value
9. **rt**: The m/z-selected rt
10. **peak_area**: Peak area of the EIC

Regardless of how the filtered spreadsheet is obtained, the path to it must be explicitly indicated within the current version of `query.py`:
```python
# Enter path to the spreadsheet containing processed and filtered spectral features
input_file = "filtered_input_file.xlsx"
```
When `query.py` is executed with the path to a valid input spreadsheet, the user will be prompted with the following list of choices:
```
Please enter the number corresponding to the desired matching criteria:
	> 1. m/z
	> 2. m/z and rt
	> 3. m/z and ccs
	> 4. m/z, rt, and ccs
	> 5. m/z, rt, ccs, and MS/MS similarity score
```
If the user selects **5**, the following options appear:
```
Please select the desired reference MS/MS database:
	> 1. Experimental
	> 2. Theoretical
```
Here, the user must decide whether to score their experimental MS/MS data against experimentally collected spectra, or against theoretical (*in silico*) spectra. The script, of course, may be run consecutively if both outputs are desired. 

### Returns
* **If the user selects options 1 - 4**: Excel spreadsheet containing potential spectral annotations, queried against the QAC database
* **If the user selects option 5**: Excel spreadsheet containing potential spectral annotations, queried against the QAC database, as well as ranked matches (if there are multiple annotations for a single spectral feature) on the basis of MS/MS cosine similarity score 
* **If the user selects option 5**: Folder ("Fragmentation Spectra") containing .png images of mirror plots (experimental vs. reference MS/MS spectra) 

### Notes
There are multiple hardset parameters within the current version of `query.py` that, if desired, can be manually overriden. The default criteria for matches within the QAC Database are: +/- 0.025 Da, +/- 0.2s, and +/- 3% for CCS, respectively. These parameters can be changed by directly manipulating the following lines of code:
```python
# Set the default parameters to +/- 0.01 Da, +/- 0.1s, and +/- 1% for CCS, respectively
c.execute('''SELECT DISTINCT compound_name FROM qacs_rt_ccs WHERE exact_mz BETWEEN ? AND ? AND rt BETWEEN ? AND ? AND average_ccs BETWEEN ? AND ? AND gradient = ? AND ccs_calibrant = ? AND column_type = ?''', (mz - 0.01, mz + 0.01, rt - 0.1, rt + 0.1, ccs * 0.99, ccs * 1.01, gradient, ccs_calibrant, column_type))
```
If the user selections option **5**, the script will automatically extract MS/MS spectra from the list of Waters .raw files. This process occurs on the basis of *either* rt and dt (if both parameters are available), or only rt (if a dt value could not be extracted from the data). The default settings for extraction are hardset at +/1 0.01s and +/1 0.1 ms around the extracted spectral feature for rt and dt, respectively. These parameters may be changed by directly manipulating the following lines of code within `query.py`:
```python
# Extract MS/MS spectrum for the experimental spectral feature with potential matches
# Set rt parameters at +/- 0.1s and dt parameters at +/- 0.5 ms
m, i = rdr.get_spectrum(1, rt - 0.1, rt + 0.1, dt_min = dt - 0.5, dt_max = dt + 0.5)
```

The current version of `query.py` also utilizes a m/z tolerance of 0.025 Da when calculating the cosine similarity score between experimental (extracted) and reference spectra. If extracted and reference m/z values are within 0.025 Da of each other, they are appended to their respective vectors and subsequently contribute to the overall similarity score between the two spectra. This threshold may be changed by directly manipulating the following line of code within `query.py`:
```python
# Change the m/z tolerance between extracted and reference to 0.01 Da
...	
	if np.abs(reference_mz[i] - extracted_mz[j]) < 0.01:
		...
```



