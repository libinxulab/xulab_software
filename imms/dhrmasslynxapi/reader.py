"""
    dhrmasslynxapi/reader.py
    Dylan H. Ross
    2018/07/24
    
    description:
        Simple API wrapper around the MassLynx python SDK for convenience. The MassLynxReader is a single class that
        wraps the MassLynxRawInfoReader, MassLynxRawChromatogramReader, and MassLynxRawScanReader objects from the 
        MassLynx SDK (version 4.6.0)
"""


from dhrmasslynxapi.sdk.MassLynxRawInfoReader import MassLynxRawInfoReader
from dhrmasslynxapi.sdk.MassLynxRawChromatogramReader import MassLynxRawChromatogramReader
from dhrmasslynxapi.sdk.MassLynxRawScanReader import MassLynxRawScanReader
from dhrmasslynxapi.sdk.MassLynxRawReader import MassLynxException


class MassLynxReader():
    """
MassLynxReader
    description:
        Simple API wrapper around Waters MassLynx python SDK
"""

    def __init__(self, raw_path):
        """
MassLynxReader.__init__
    description:
        Initializes a new MassLynxReader from a path to the raw data file
    parameters:
        raw_path (str) -- path to the raw data file
"""
        # store the path to the raw file and initialize a reader object
        self.raw = raw_path
        
        # initialize info, chromatogram, and scan reader objects from the SDK
        self.info_reader = MassLynxRawInfoReader(raw_path)
        self.chrom_reader = MassLynxRawChromatogramReader(self.info_reader)
        self.scan_reader = MassLynxRawScanReader(self.info_reader)

        # establish some basic information about the data file
        self.n_funcs = self.info_reader.GetNumberofFunctions()
        self.scans_per_func = [self.info_reader.GetScansInFunction(_) for _ in range(self.n_funcs)]
        time_ranges_per_func = [self.info_reader.GetAcquisitionTimeRange(_) for _ in range(self.n_funcs)]
        # arrays mapping scan numbers (indices) to scan times (values)
        self.scan_times = [[i * (r[1] - r[0]) / (n - 1) + r[0] for i in range(n)] for n, r in zip(self.scans_per_func, 
                                                                                                  time_ranges_per_func)]
            

    def __repr__(self):
        """
MassLynxReader.__repr__
    description:
        returns a string representation of this MassLynxReader, consists of 
        object name and path to raw file
    returns:
        (str) -- string representation
"""
        return "MassLynxReader({})".format(self.raw)

  
    def get_chrom(self, func, mass, tol):
        """
MassLynxReader.get_chrom
    description:
        Returns a chromatogram from the specified function (e.g. LC, IM) corresponding to a target mass and tolerance. 
        Returns separate lists of times and intensities, or None, None if anything goes wrong
    parameters:
        func (int) -- function to use for extracting chromatogram
        mass (float) -- target mass
        tol (float) -- mass tolerance
    returns:
        tuple(list(float), list(float)) -- times, intensities
            OR
        tuple(None, None) -- if an error occurs
"""
        try:
            times, intensities = self.chrom_reader.ReadMassChromatogram(func, mass, tol, False)
            return (times, intensities)
        except Exception as e:
            print("MassLynxReader: get_chrom: failed to generate chromatogram ({})\n    file: {}\n    function: {}" \
                "\n    mass: {:.4f}".format(e, self.raw, func, mass))
            return (None, None)


    def __func_has_mobility(self, func):
        """ checks whether a given function contains mobility data """
        has_mobility = False
        try:
            self.info_reader.GetDriftTime(func, 0)  # will throw an error if not drift time in this function
            has_mobility = True
        except MassLynxException:
            pass
        return has_mobility


    def __get_scan_indices(self, func, t_min, t_max):
        """ uses min and max scan times to compute corresponding scan indices """
        scan_indices = []
        for i in range(self.scans_per_func[func]):
            t_i = self.scan_times[func][i] 
            if t_i >= t_min and t_i <= t_max:
                scan_indices.append(i)
        if scan_indices == []:
            # if we did not get scan indices then the window may have been too small
            m = 'MassLynxReader: __get_scan_indices: unable to compute scan indices for bounds: ' \
                '({}, {})'.format(t_min, t_max)
            raise ValueError(m)
        return scan_indices


    def __get_drift_indices(self, func, dt_min, dt_max):
        """ uses min and max drift times to compute corresponding drift indices """
        drift_indices = []
        for i in range(200):  # can safely assume there are 200 dt bins when mobility data is available
            dt_i = self.info_reader.GetDriftTime(func, i) 
            if dt_i >= dt_min and dt_i <= dt_max:
                drift_indices.append(i)
        if drift_indices == []:
            # if we did not get drift indices then the window may have been too small
            m = 'MassLynxReader: __get_drift_indices: unable to compute drift indices for bounds: ' \
                '({}, {})'.format(dt_min, dt_max)
            raise ValueError(m)
        return drift_indices


    def get_spectrum(self, func, t_min, t_max, dt_min=None, dt_max=None, accum_precis=4):
        """
MassLynxReader.get_spectrum
    description:
        Accumulates MS scans from a given function between a set of time bounds (retention time or drift time depending
        on the function used). Returns lists of m/z values and intensities (summed from all scans) or None, None if 
        anything goes wrong
    parameters:
        func (int) -- function to use for getting the spectrum
        t_min (float) -- lower time bound (retention time or drift time)
        t_max (float) -- upper time bound (retention time or drift time)
        [dt_min, dt_max (None or float)] -- if not None, use the specified drift time bounds for filtering as well 
                                            [optional, default=None]
        [accum_precis (int)] -- mass precision (decimal places) to use when accumulating scans (a.k.a. mass binning), 
                                empirical testing with MassLynx showed 4 to be the best at recapitulating the exact
                                masses that one might observe using Spectrum  [optional, default=4]
    returns:
        tuple(list(float), list(float)) -- m/z, intensities (summed from all scans in window)
            OR
        tuple(None, None) -- if an error occurs
"""
        # determine the scan indices to use
        scan_indices = self.__get_scan_indices(func, t_min, t_max)
        
        drift_indices = []
        # check that mobility data is included in the specified function if drift time bounds are provided
        if dt_min is not None and dt_max is not None:
            if not self.__func_has_mobility(func):
                # raise exception
                m = 'MassLynxReader: get_spectrum: function {} does not appear to contain drift time'.format(func)
                raise ValueError(m)
            # compute the drift indices
            drift_indices = self.__get_drift_indices(func, dt_min, dt_max)
        elif dt_min is not None or dt_max is not None:
            # if you are going to provide one bound, you need to provide both
            m = 'MassLynxReader: get_spectrum: both drift time bounds must be provided to filter on drift time'
            raise ValueError(m)

        # get the scan data
        if drift_indices != []:
            scans = [self.scan_reader.ReadDriftScan(func, s, d) for s in scan_indices for d in drift_indices]
        else:
            scans = [self.scan_reader.ReadScan(func, s) for s in scan_indices]

        # accumulate scans into single spectrum
        d = {}
        for scan in scans:
            for n in range(len(scan[0])):
                m = round(scan[0][n], accum_precis)
                i = scan[1][n]
                if m in d:
                    d[m] += i
                else:
                    d[m] = i
        m = []
        i = []
        for m_, i_ in sorted(d.items()):
            m.append(m_)
            i.append(i_)

        return m, i


    def get_filtered_chrom(self, func, mass, tol, rt_min=None, rt_max=None, dt_min=None, dt_max=None):
        """
MassLynxReader.get_filtered_chrom
    description:
        Returns a chromatogram from the specified function corresponding to a target mass and tolerance, but with
        filtering based on retention time or drift time. The specified function must contain mobility data.

        * This method is much slower than the get_chrom method, so it is preferred to use it only when you need to 
          filter data based on retention time or drift time *

        * Only one set of bounds (retention time or drift time) can be used at once *
    parameters:
        func (int) -- function to use for generating chromatogram
        mass (float) -- target mass
        tol (float) -- mass tolerance
        [rt_min, rt_max (None or float)] -- if not None, use the specified retention time bounds for filtering 
                                            [optional, default=None]
        [dt_min, dt_max (None or float)] -- if not None, use the specified drift time bounds for filtering
                                            [optional, default=None]
    returns:
        tuple(list(float), list(float)) -- times, intensities
            OR
        tuple(None, None) -- if an error occurs
"""
        # make sure that the function has mobility data
        if not self.__func_has_mobility(func):
            m = 'MassLynxReader: get_filtered_chrom: function {} does not appear to contain drift time'.format(func)
            raise ValueError(m)

        # make sure only one set of bounds is provided
        rt_bounds_provided = rt_min is not None and rt_max is not None
        dt_bounds_provided = dt_min is not None and dt_max is not None
        if rt_bounds_provided and dt_bounds_provided:
            m = 'MassLynxReader: get_filtered_chrom: only rt or dt bounds can be provided at once' 
            raise ValueError(m)
        # not enough bounds provided 
        if not rt_bounds_provided and not dt_bounds_provided:
            m = 'MassLynxReader: get_filtered_chrom: complete rt or dt bounds need to be provided'
            raise ValueError(m)

        scans, times, intensities = [], [], []

        # retention time chromatogram
        if rt_bounds_provided:
            # maximal dt bounds
            drift_indices = [_ for _ in range(0, 200)]
            scan_indices = self.__get_scan_indices(func, rt_min, rt_max)
            scans = [[self.scan_reader.ReadDriftScan(func, s, d) for s in scan_indices] for d in drift_indices]
            # times is drift time
            times = [self.info_reader.GetDriftTime(func, i) for i in drift_indices]

        # drift time chromatogram
        if dt_bounds_provided:
            drift_indices = self.__get_drift_indices(func, dt_min, dt_max)
            # maximal rt bounds
            scan_indices = [_ for _ in range(self.scans_per_func[func])]
            scans = [[self.scan_reader.ReadDriftScan(func, s, d) for d in drift_indices] for s in scan_indices]
            # times is retention time
            times = self.scan_times[func]

        for dim1_scans in scans:
            acc_int = 0.  # accumulated intensity over scans
            for scan in dim1_scans:
                for mz, inten in zip(*scan):
                    if abs(mass - mz) <= tol:
                        acc_int += inten
            intensities.append(acc_int)

        return times, intensities
        


