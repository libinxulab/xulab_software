"""
    dhrmasslynxapi/ccs_calibration.py
    Dylan Ross
    2018/08/09

    description:
        CCSCalibrationBase is a convenient object for performing CCS calibration
        akin to the core module found in CcsCal. I am opting for a lighter,
        purpose-built reimplementation rather than shoehorning in the
        implementation from CcsCal.
"""


import matplotlib
matplotlib.use('AGG')
from scipy.optimize import curve_fit
from numpy import exp, sqrt, argmax, linspace
from matplotlib import pyplot as plt
from matplotlib import gridspec as gs
from pandas import read_excel, ExcelWriter
import os


from dhrmasslynxapi.reader import MassLynxReader


class CCSCalibrationBase:
    """
CCSCalibrationBase
    description:
        Modular base object for making and applying CCS calibrations
"""

    def __init__(self):
        """
CCSCalibrationBase.__init__
    description:
        TODO
"""
        # individual calibrant information stored here
        self.calibrants = []
        # calibration curve parameters
        self.cal_params = (None, None, None)
        # masses, dts, and ccss used to fit the calibration curve
        self.cal_data = (None, None, None)
        self.cal_corr_data = (None, None)

    @staticmethod
    def reduced_mass(mass, ref_mass=28.0134):
        """
CCSCalibrationBase.reduced_mass
    description:
        Calculates reduced mass of an ion using a reference mass
    parameters:
        mass (float) -- input mass
        [ref_mass (float)] -- reference mass, defaults to N2
                                [optional, default=28.0134]
    returns:
        (float) -- reduced mass
"""
        return (mass * ref_mass) / (mass + ref_mass)

    @staticmethod
    def corrected_dt(dt, mass, edc=1.35):
        """
CCSCalibrationBase.corrected_dt
    description:
        Calculates a drift time corrected for mass-dependent flight time
    parameters:
        dt (float) -- original drift time
        mass (float) -- the m/z to use for the correction
        [edc (float)] -- EDC delay coefficient [optional, default=1.35]
    returns:
        (float) -- corrected drift time
"""
        return dt - (sqrt(mass) * edc) / 1000.

    @staticmethod
    def corrected_ccs(ccs, mass, charge):
        """
CCSCalibrationBase.corrected_ccs
    description:
        Calculates CCS corrected for mass-dependent flight time
    parameters:
        ccs (float) -- original CCS
        mass (float) -- the m/z to use for the correction
        charge (float) -- charge state
    returns:
        (float) -- corrected CCS
"""
        return ccs * (charge * sqrt(CCSCalibrationBase.reduced_mass(mass)))

    @staticmethod
    def cal_curve(dt_c, A, t0, B):
        """
CCSCalibrationBase.cal_curve
    description:
        Basic power function for calibration curve of the form:
            CCS' = A(dt' + t0)^B
    parameters:
        dt_c (float) -- corrected drift time
        A (float) -- A parameter
        t0 (float) -- t0 parameter
        B (float) -- B parameter
    Returns:
        (float) -- CCS' (divide by reduced mass for actual CCS...)
    """
        return A * (dt_c + t0)**B

    def fit_cal_curve(self, masses, dts, ccss):
        """
CCSCalibrationBase.fit_cal_curve
    description:
        Fits a calibration curve to dt and CCS values and stores the optimized
        parameters in self.cal_params. Expects UNCORRECTED dts and ccss
    parameters:
        masses (list(float)) -- calibrant masses
        dts (list(float)) -- uncorrected calibrant drift times
        ccss (list(float)) -- uncorrected calibrant ccss
"""
        corr_dt = [self.corrected_dt(dt, mass) for dt, mass in zip(dts, masses)]
        corr_ccs = [self.corrected_ccs(ccs, mass, self.charge) for ccs, mass in zip(ccss, masses)]
        self.cal_params, cov = curve_fit(self.cal_curve,
                                         corr_dt, corr_ccs,
                                         maxfev=1000000, p0=(500., 0.001, 0.5))
        # after a successful fit store the masses, drift times, ccss used to
        # create the calibration curve
        self.cal_data = (masses, dts, ccss)
        self.cal_corr_data = (corr_dt, corr_ccs)

    def calibrated_ccs(self, mass, dt):
        """
CCSCalibrationBase.calibrated_ccs
    description:
        Use the calibration curve parameters in self.cal_params to
        compute calibrated CCS for an m/z dt pair
    parameters:
        mass (float) -- m/z
        dt (float) -- drift time
    returns:
        (float) -- calibrated ccs
"""
        if not self.cal_params[0]:
            raise RuntimeError("CCSCalibrationBase: calibrated_ccs: self.cal_params not set")
        return self.cal_curve(self.corrected_dt(dt, mass), *self.cal_params) / (self.charge * sqrt(self.reduced_mass(mass)))

    def cal_curve_figure(self, fig_name):
        """
CCSCalibrationBase.cal_curve_figure
    description:
        Produces a figure from the CCS calibration curve with residuals and
        saves it as a .png
    parameters:
        fig_name (str) -- filename to save the calibration curve under
"""
        if not self.cal_params[0]:
            raise RuntimeError("CCSCalibrationBase: cal_curve_figure: self.cal_params not set")
        # gather necessary data
        corr_dt, corr_ccs = self.cal_corr_data
        masses, dts, ccss = self.cal_data
        corr_ccs_calc = [self.corrected_ccs(self.calibrated_ccs(mass, dt), mass, self.charge) for mass, dt in zip(masses, dts)]
        ccs_resid = [100. * (ccs - self.calibrated_ccs(mass, dt)) / ccs for mass, dt, ccs in zip(masses, dts, ccss)]
        # set up the figure
        fig = plt.figure(figsize=(5, 4.5))
        grid = gs.GridSpec(2, 1, height_ratios=(5, 2))
        ax1, ax2 = fig.add_subplot(grid[0]), fig.add_subplot(grid[1])
        ax2.axhline(lw=0.75, ls="--", c="k")
        # plot the data
        ax1.plot(corr_dt, corr_ccs, "bo", ms=4, mew=1, fillstyle='none', label='calibrants')
        x = linspace(min(corr_dt), max(corr_dt), 100)
        y = self.cal_curve(x, *self.cal_params)
        ax1.plot(x, y, "b-", lw=1, alpha=0.6, label='fitted')
        ax1.legend(frameon=False)
        ax2.bar(corr_dt, ccs_resid, 0.25, color=(0., 0., 1., 0.6), align='center')
        # axis adjustments and labels
        ax1.set_title("CCS calibration")
        ax1.set_ylabel("corrected CCS")
        ax2.set_ylabel("residual CCS (%)")
        ax2.set_xlabel("corrected drift time (ms)")
        plt.savefig(fig_name, dpi=400, bbox_inches='tight')
        plt.close()

    def __str__(self):
        """
CCSCalibrationBase.__str__
    description:
        produces a string representation of this instance, contains information
        about the calibration curve (if fitted)
    returns:
        (str) -- string representation of this CCSCalibrationBase object
"""
        if not self.cal_params[0]:
            return "calibration curve not fitted"
        # calibration has been completed, include the fitted calibration curve parameters
        else:
            return "fitted parameters: \n\tA = {:.3f}\n\tt0 = {:.3f}\n\tB = {:.3f}".format(*self.cal_params)


class CCSCalibrationList(CCSCalibrationBase):
    """
CCSCalibrationList
    description:
        A convenient object for making and applying CCS calibrations from lists of m/z, dt, and ref CCS
"""

    def __init__(self, masses, dts, ccss, charge=1.):
        """
CCSCalibrationList.__init__
    description:
        Initializes a new CCSCalibrationList object and constructs a calibration
        curve based on lists of m/z, dt, and ref CCS
    parameters:
        masses (list(float)) -- m/z for calibrants
        dts (list(float)) -- drift times for calibrants
        ccss (list(float)) -- reference CCS for calibrants
        [charge (float)] -- charge state [optional, default=1.] 
"""
        # call the superclass initializer
        super().__init__()
        
        # store charge
        self.charge = charge
        
        # store calibrant data
        for m, d, c in zip(masses, dts, ccss):
            self.calibrants.append({
                'mz': m,
                'dt': d,
                'ref_ccs': c,
            })

        # construct the CCS calibration
        self.fit_cal_curve(masses, dts, ccss)

        # add residual ccs (%) attributes to all of the calibrants
        for calibrant in self.calibrants:
            calibrant['resid_ccs'] = 100. * (calibrant['ref_ccs'] - self.calibrated_ccs(calibrant['mz'], calibrant['dt'])) / calibrant['ref_ccs']


class CCSCalibrationRaw(CCSCalibrationList):
    """
CCSCalibrationRaw
    description:
        A convenient object for making and applying CCS calibrations directly
        from raw data files.
"""

    def __init__(self, raw, masses, ccss, mass_window=0.05, dt_func=1, no_init=False, charge=1.):
        """
CCSCalibrationRaw.__init__
    description:
        Initializes a new CCSCalibrationRaw object and constructs a calibration
        curve based on calibrants with specified masses from the specified data
        file(s)
    parameters:
        raw (str OR
             list(str)) -- single data file (str) or multiple data files
                            (list(str)) to extract data from
        masses (list(float) OR
                list(list(float))) -- single list (list(float)) of masses or
                                        multiple lists (list(list(float)))
                                        of masses for calibrants
        ccss (list(float) OR
              list(list(float))) -- single list (list(float)) of reference CCS
                                    values or multiple lists (list(list(float)))
                                    of reference CCS values for calibrants
        [mass_window (float)] -- mass tolerance [optional, default=0.05]
        [dt_func (int)] -- the MS function containing drift time data
                            [optional, default=1]
        [no_init (bool)] -- do not perform any initialization, just make a blank CCSCalibrationRaw instance
                            [optional, default=False]
        [charge (float)] -- charge state [optional, default=1.] 
"""
        # make sure no_init was not set
        if not no_init:
            if type(raw) == str:
                # for single raw file and mass/ccs list, convert to lists with those single values then proceed
                raw, masses, ccss = [raw], [masses], [ccss]
            # iterate through the raw file and mass list combinations
            # and store detailed information on each calibrant
            raws, dts, atds, gauss_params, masses2, ccss2 = [], [], [], [], [], []
            for raw_, masses_, ccss_ in zip(raw, masses, ccss):
                print("processing {} for masses: {}".format(raw_, masses_))
                reader = MassLynxReader(raw_)
                for mass, ccs in zip(masses_, ccss_):
                    # extract the ATD
                    drift_time, intensity = reader.get_chrom(dt_func, mass, mass_window)
                    # raise an error if we did not get data out for some reason
                    if not drift_time:
                        e_msg = "CCSCalibrationRaw: __init__: unable to fit drift time"
                        raise RuntimeError(e_msg)
                    # fit the ATD with a gaussian to get the drift time (parameter B)
                    A, B, C = self.peak_fit(drift_time, intensity)
                    print("\tmz {:.4f} dt {:.2f}".format(mass, B))
                    # record calibrant data and metadata
                    dts.append(B)
                    raws.append(raw)
                    atds.append((drift_time, intensity))
                    gauss_params.append((A, B, C))
                    masses2.append(mass)
                    ccss2.append(mass)

            # create calibration using superclass __init__
            super().__init__(masses2, dts, ccss2, charge=charge)

            # add calibrant metadata from drift time extraction to self.calibrants
            for calibrant, raw, atd, gp in zip(self.calibrants, raws, atds, gauss_params):
                calibrant['raw'] = raw
                calibrant['atd'] = atd
                calibrant['gauss_params'] = gp


    @staticmethod
    def gauss(x, A, B, C):
        """
CCSCalibrationRaw.gauss
    description:
        basic gaussian function for peak fitting
    parameters:
        x (float) -- input
        A (float) -- peak height (amplitude)
        B (float) -- peak center (mean)
        C (float) -- peak width (standard deviation)
    returns:
        (float) -- output
"""
        if abs(C) < 0.01:
            C = 0.01 
        return A * exp(-(x - B)**2 / (2. * C**2))

    def peak_fit(self, drift_time, intensity, p0='guess'):
        """
CCSCalibrationRaw.peak_fit
    description:
        Fits a Gaussian function to a supplied ATD, returning the fitted
        parameters.
    parameters:
        drift_time (list(float)) -- drift time component of the ATD
        intensity (list(float)) -- intensity component of the ATD
        [p0 (str or tuple(float))] -- initial parameters for peak fit may be guessed ('guess') or provided directly
                                        [optional, default='guess']
    returns:
        tuple(float) -- tuple containing fitted Gaussian A, B, C parameters
"""
        if p0 == 'guess':
            p0 = (max(intensity), drift_time[argmax(intensity)], 0.5)
        opt, cov = curve_fit(self.gauss,
                             drift_time, intensity,
                             maxfev=5000,
                             p0=p0)
        return opt

    def __str__(self):
        """
CCSCalibrationRaw.__str__
    description:
        produces a string representation of this instance, contains information
        about the calibration curve (if fitted) and the raw file(s) used to
        generate the calibration
    returns:
        (str) -- string representation of this CCSCalibrationRaw object
"""
        # make raw file list
        repr = "CCS calibration from RAW files:\n"
        for raw in self.raw:
            repr += "\t" + raw + "\n"
        # add in CCS calibration curve parameters from the superclass and return
        return repr + super().__str__()


class CCSCalibrationRawXL(CCSCalibrationRaw):
    """
CCSCalibrationRawXL
    description:
        A convenient object for making and applying CCS calibrations directly
        from raw data files that takes a single .xlsx spreadsheet as input where
        the sheet names correspond to the data files and each sheet consists of
        a list of calibrant masses and reference CCS values. Requires pandas and
        xlrd libraries.
"""

    def __init__(self, xlsx, mass_window=0.05, dt_func=1, charge=1., path_prefix="", name_prefix=""):
        """
CCSCalibrationRawXL.__init__
    description:
        Initializes a new CCSCalibrationRaw object and constructs a calibration
        curve based on calibrants with specified masses from the specified data
        file(s), all specified in a single Excel .xlsx spreadsheet
    parameters:
        xlsx (str) -- Excel .xlsx spreadsheet to read calibrant info from
        [mass_window (float)] -- mass tolerance [optional, default=0.05]
        [dt_func (int)] -- the MS function containing drift time data
                            [optional, default=1]
        [charge (float)] -- charge state [optional, default=1.] 
        [path_prefix (str)] -- prefix to prepend to paths to data files, default behavior is to
                                look in the same directory as the input xlsx file if no path_prefix
                                is provided [optional, default=""]
        [name_prefix (str)] -- prefix to prepend to raw file names taken from the xlsx input
                                sheet names [optional, default=""]
"""
        # load the excel file
        xl = read_excel(xlsx, sheet_name=None, index_col=None, header=None)
        # get the path to the raw files (same directory as input xlsx, OR specified in path_prefix)
        if not path_prefix:
            raw_dir = os.path.split(xlsx)[0]
        else:
            raw_dir = path_prefix
        raws = [os.path.join(raw_dir, name_prefix + sheet) for sheet in xl]
        masses = [xl[sheet].values.T.tolist()[0] for sheet in xl]
        ccss = [xl[sheet].values.T.tolist()[1] for sheet in xl]
        # call the superclass initializer
        super().__init__(raws, masses, ccss, mass_window=mass_window, dt_func=dt_func, charge=charge)

    def batch_calibrated_ccs(self, xlsx, xlsx_out=None):
        """
CCSCalibrationRawXL.batch_calibrated_ccs
    description:
        get calibrated CCS for a set of m/z and drift time values provided in
        an Excel .xlsx spreadsheet. The Excel spreadsheet must contain columns
        with the headers "file", "m/z" and "drift time". A column with the header 'ccs'
        will be written to the spreadsheet, and if one already exists it will
        be overridden. Requires openpyxl for ExcelWriter.
    parameters:
        xlsx (str) -- Excel .xlsx spreadsheet to read input m/z and drift times
                        from
        [xlsx_out (str)] -- file name to save the output .xlsx file under
                            [optional, default=xlsx]
"""
        if not self.cal_params[0]:
            raise RuntimeError("CCSCalibrationRawXL: batch_calibrated_ccs: self.cal_params not set")
        # load the excel file
        xl = read_excel(xlsx, sheet_name=None, index_col=None)
        # set the output file name, defaults to same as input
        if not xlsx_out:
            xlsx_out = xlsx
        # output writer
        xl_w = ExcelWriter(xlsx_out)
        # process sheet by sheet
        for sheet in xl:
            # if a sheet does not contain 'm/z' or 'drift time'...
            # ... then skip it? or throw an error?
            for column in ('m/z', 'drift time'):
                if column not in xl[sheet]:
                    msg = "CCSCalibrationRawXL: batch_calibrated_ccs: column {} not found in sheet '{}'".format(column, sheet)
                    # throw an error because I am unforgiving
                    raise ValueError(msg)
            masses = xl[sheet]['m/z'].tolist()
            dts = xl[sheet]['drift time'].tolist()
            print(masses)
            print(dts)
            ccss = [self.calibrated_ccs(mass, dt) for mass, dt in zip(masses, dts)]
            # add the ccs column to the sheet
            xl[sheet]['ccs'] = ccss
            # write the sheet to the output file
            xl[sheet].to_excel(xl_w, sheet_name=sheet, index=False)
        # close the ExcelWriter
        xl_w.save()


