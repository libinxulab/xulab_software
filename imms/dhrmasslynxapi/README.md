## `dhrmasslynxapi`

This Python package is a wrapper around the MassLynx API from Waters, and is used for extracting data from 
Waters .raw files. There is also a module for performing CCS calibration. This Package requires the original Waters
sdk files to be copied into the `sdk/` directory, which are not included in this repository but can be obtained from
the lab NAS (under `lab_resources/masslynx_sdk_files/`):
* `cdt.dll`
* `MassLynxLockMassProcessor.py`
* `MassLynxParameters.py`
* `MassLynxRaw.dll`
* `MassLynxRaw.lib`
* `MassLynxRawChromatogramReader.py`
* `MassLynxRawDefs.py`
* `MassLynxRawInfoReader.py`
* `MassLynxRawReader.py`
* `MassLynxRawScanReader.py`
* `MassLynxSampleList.py`


## `dhrmasslynxapi.reader`
This module defines the `MassLynxReader` object which is the primary object for interacting with .raw MS data files. 


### Initialization
The `MassLynxReader` object is initialized with the path to the .raw file:
```python
from dhrmasslynxapi.reader import MassLynxReader

rdr = MassLynxReader('path/to/my/data/file.raw')
```
Upon initializing an instance of `MassLynxReader`, instance variables are created describing the functions contained in
the data file:
* `MassLynxReader.n_funcs` - number of functions
* `MassLynxReader.scans_per_func` - number of scans in each function
* `MassLynxReader.scan_times` - arrays mapping scan numbers (indices) to scan times (values)


### Methods

#### `MassLynxReader.get_chrom`
Returns a chromatogram from the specified function (_e.g._ LC, IM) corresponding to a target mass and tolerance. 
Returns separate lists of times and intensities, or `(None, None)` if anything goes wrong

##### Parameters
* `func` (`int`) - function to use for extracting chromatogram
* `mass` (`float`) - target m/z
* `tol` (`float`) - mass tolerance

##### Returns
* `tuple(list(float), list(float))` - lists of times and intensities, _i.e._ the extracted chromatogram


#### `MassLynxReader.get_filtered_chrom`
Returns a chromatogram from the specified function (_e.g._ LC, IM) corresponding to a target mass and tolerance, 
filtering based on retention time or drift time. The specified function must contain mobility data.
  
  
*This method is much slower than the `get_chrom` method, so it is preferred to use it only when you need to filter 
data based on retention time or drift time*.
  

*Only one set of bounds (retention time or drift time) can be used at once*. 
  

Returns separate lists of times and intensities, or `(None, None)` if anything goes wrong

##### Parameters
* `func` (`int`) - function to use for extracting chromatogram
* `mass` (`float`) - target m/z
* `tol` (`float`) - mass tolerance
* [`rt_min`, `rt_max` (`None` or `float`)] - if not `None`, use the specified retention time bounds for 
filtering [optional, default=`None`]
* [`dt_min`, `dt_max` (`None` or `float`)] - if not `None`, use the specified drift time bounds for 
filtering [optional, default=`None`]

##### Returns
* `tuple(list(float), list(float))` - lists of times and intensities, _i.e._ the extracted chromatogram

##### Example
```python
# extract LC chromatogram from function 1 (includes mobility data) for m/z 234.5678 +/- 0.01 
# and having drift time between 3.4 and 4.5 ms
rt, i = rdr.get_filtered_chrom(1, 234.5678, 0.01, dt_min=3.4, dt_max=4.5)
```


#### `MassLynxReader.get_spectrum`
Accumulates MS scans from a given function between a set of time bounds (retention time or drift time depending on the 
function used). Returns lists of m/z values and intensities (summed from all scans) or `(None, None)` if  anything 
goes wrong

##### Parameters
* `func` (`int`) - function to use for extracting the spectrum
* `t_min` (`float`) -- lower time bound (retention time or drift time, depending on the function)
* `t_max` (`float`) -- upper time bound (retention time or drift time, depending on the function)
* [`dt_min`, `dt_max` (`None` or `float`)] -- if not `None`, use the specified drift time bounds for filtering in 
addition to retention time, requires mobility data in the specified function [optional, default=`None`]
* [`accum_precis` (`int`)] -- mass precision (decimal places) to use when accumulating scans (_a.k.a._ mass binning), 
empirical testing with MassLynx showed `4` to be the best at recapitulating the exact masses that one might observe 
using the Spectrum program [optional, default=`4`]

##### Returns
* `tuple(list(float), list(float))` - m/z and intensities (summed from all scans in window), _i.e._ the 
extracted spectrum

##### Example
```python
# extract mass spectrum from function 1 (includes mobility data) for retention time range 1.2 to 2.3 min
# and drift time range 3.4 to 4.5 ms
m, i = rdr.get_spectrum(1, 1.2, 2.3, dt_min=3.4, dt_max=4.5)
```


## `dhrmasslynxapi.ccs_calibration`
This module defines three objects `CCSCalibrationList`, `CCSCalibrationRaw` and `CCSCalibrationRAWXL`. Upon 
initialization, these objects construct a CCS calibration curve (using calibrant data provided either in the form of 
lists or extracted automatically from raw data files). After initialization, calibrated CCS values can be obtained for
mass and drift time combinations using the `calibrated_ccs` method. 


### Initialization
The basic CCS calibration object (`CCSCalibrationList`) can be initialized using lists of masses, extracted drift 
times, and reference CCS values for all calibrants. 
CCS calibration objects can also be initialized using either the path to .raw data file(s) and a list of masses and 
reference CCS values for CCS calibrants (`CCSCalibrationRaw`) or using the path to an excel spreadsheet containing the 
same information (`CCSCalibrationRawXL`). Both CCS calibration objects can be set to extract calibrant data from one 
data file or multiple data files. When the CCS calibration object is initialized, it automatically extracts the drift
times from calibrant data file(s), then fits a CCS calibration curve using the reference m/z and CCS values and
these fitted drift times. 


##### Parameters (`CCSCalibrationList`)
* `masses` (`list(float)`) - m/z for calibrants
* `dts` (`list(float)`) - drift times for calibrants
* `ccss` (`list(float)`) - reference CCS for calibrants
* [`charge` (`float`)] - charge state [optional, default=`1.`] 

##### Example (`CCSCalibrationList`)
```python
from dhrmasslynxapi.ccs_calibration import CCSCalibrationList

# m/z, drift time, and reference CCS values for polyalanine CCS calibrants
polyala_mz = [123.4567, 234.5678, ...]
polyala_dt = [12.3, 23.4, ...]
polyala_ccs = [123.4567, 234.5678, ...]

cal_polyala = CCSCalibrationList(polyala_mz, polyala_dt, polyala_ccs)

```


##### Parameters (`CCSCalibrationRaw`)
* `raw` (`str` or `list(str)`) - single data file (str) or multiple data files(list(str)) to extract data from
* `masses` (`list(float)` or `list(list(float))`) - single list (`list(float)`) of masses or multiple lists 
(`list(list(float))`) or masses for calibrants
* `ccss` (`list(float)` or `list(list(float))`) - single list (`list(float)`) of reference CCS values or multiple 
lists (`list(list(float))`) of reference CCS values for calibrants
* [`mass_window` (`float`)] - mass tolerance [optional, default=`0.05`]
* [`dt_func` (`int`)] - the MS function containing drift time data [optional, default=`1`]
* [`no_init` (`bool`)] - do not perform any initialization, just make a blank CCSCalibrationRaw instance [optional, 
default=`False`]
* [`charge` (`float`)] - charge state [optional, default=`1.`] 

##### Example (`CCSCalibrationRaw`)
```python
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRaw

# m/z and reference CCS values for polyalanine CCS calibrants
polyala_mz = [123.4567, 234.5678, ...]
polyala_ccs = [123.4567, 234.5678, ...]

# m/z and reference CCS values for drug CCS calibrants
drug_mz = [123.4567, 234.5678, ...]
drug_ccs = [123.4567, 234.5678, ...]

# initialize CCS calibration object (polyalanine CCS calibrants)
cal_polyala = CCSCalibrationRaw('polyala.raw', polyala_mz, polyala_ccs)

# initialize CCS calibration object (drug CCS calibrants)
cal_drug = CCSCalibrationRaw('drugs.raw', drug_mz, drug_ccs)

# initialize CCS calibration object (polyalanine and drug CCS calibrants)
cal_polyala_drug = CCSCalibrationRaw(['polyala.raw', 'drugs.raw'], 
                                     [polyala_mz, drug_mz], 
                                     [polyala_ccs, drug_ccs])

# initialize CCS calibration object (polyalanine CCS calibrants, combine data from 3 replicates)
cal_polyala_3reps = CCSCalibrationRaw(['polyala_rep1.raw', 'polyala_rep2.raw', 'polyala_rep3.raw'], 
                                      [polyala_mz for _ in range(3)], 
                                      [polyala_ccs for _ in range(3)])
```


##### Parameters (`CCSCalibrationRawXL`)
* `xlsx` (`str`) -- Excel .xlsx spreadsheet to read calibrant info from. Sheet names correspond to the calibrant data
files, and each sheet contains calibrant m/z and reference CCS values
* [`mass_window` (`float`)] - mass tolerance [optional, default=`0.05`]
* [`dt_func` (`int`)] - the MS function containing drift time data [optional, default=`1`]
* [`no_init` (`bool`)] - do not perform any initialization, just make a blank CCSCalibrationRaw instance [optional, 
default=`False`]
* [`charge` (`float`)] - charge state [optional, default=`1.`]
* [`path_prefix` (`str`)] - prefix to prepend to paths to data files, default behavior is to look in the same 
directory as the input xlsx file if no `path_prefix` is provided [optional, default=`""`]
* [`name_prefix` (`str`)] - prefix to prepend to raw file names taken from the xlsx input sheet 
names [optional, default=`""`]

##### Example (`CCSCalibrationRawXL`)
```python
from dhrmasslynxapi.ccs_calibration import CCSCalibrationRawXL

# initialize CCS calibration object (polyalanine CCS calibrants, combine data from 3 replicates)
cal_polyala_3reps = CCSCalibrationRawXL('cal_polyala_3reps.xlsx')
```

_In this example`cal_polyala_3reps.xlsx` is an excel spreadsheet where each sheet contains the masses and literature 
CCS values for a set of calibrants to extract from a single data file, which is specified by the sheet name. To extract 
calibrant data from multiple data files, add corresponding sheets for each data file that contain the respective 
calibrant information. The first column in each sheet is the calibrant m/z and the second is the literature CCS, there
should be no headers_


### `calibrated_ccs`
This is the primary method to use after initializing a CCS calibration object. Takes m/z and drift time as input and
returns a calibrated CCS value based on the calibration curve established during initialization. This is the same for
`CCSCalibrationList`, `CCSCaibrationRaw`, and `CCSCalibrationRawXL` subclasses.

##### Parameters
* `mass` (`float`) - m/z
* `dt` (`float`) - drift time

#### Returns
* `ccs` (`float`) - calibrated ccs

##### Example
```python
# cal is an initialized instance of CCSCalibrationList, CCSCalibrationRaw, or CCSCalibrationRawXL

# get a single calibrated CCS value
ccs = cal.calibrated_ccs(234.5678, 3.45)

# get calibrated CCS values for a list of unknown m/zs and drift times
unknown_mzs = [123.4567, 234.5678, ...]
unknown_dts = [1.23, 2.34, ...]
unknown_ccss = [cal.calibrated_ccs(mz, dt) for mz, dt in zip(unknown_mzs, unknown_dts)]

```
