"""serfling.py creates a serfling regression from values from an excel file.

Usage:

Running this script will estimate serfling regression parameters and other
useful data, such as std dev of the error between the data and the curve,
and other curve offsets.

> python serfling.py <file_path> <cycles>

file_path: the path of the excel file (eg. data.xlsx)
cycles: the number of cycles that the provided data has.
    used to manually calculate the period to avoid overfitting
"""

import numbers
import sys

import numpy as np

from scipy.optimize import leastsq
from scipy.stats import linregress

from openpyxl import load_workbook

def GetData(file_name, sheet=None):
    """Extracts data from the file_name xlsx file.

    Gets all the datapoints located in the A and B columns of an excel sheet,
    up until the first empty A or B cell, where A is time, and B is the value.

    Ignores the first row (Since it usually has column names), and only extracts
    values if the respective C column cell contains the string "Baseline"

    Params:
        file_name: The name of a file_name file. Should be located in the script's
            folder
        sheet: (string) name of the excel sheet from which to extract the data.
            defaults to the first sheet in the file.

    Returns:
        Two numpy arrays, where the first corresponds to the time entry for each data
        point, and the second contains the respective data values for each time
        entry.
        The lists are sorted by time.

    Raises:
        ValueError: If the provided sheet doesn't exist in the provided file
            or if a input data point is not conformed by numbers.
    """
    excel_workbook = load_workbook(file_name)
    sheets = excel_workbook.get_sheet_names()
    if sheet:
        if sheet not in sheets:
            print "ERROR: Sheet", sheet, "does not exist in file", file_name
            raise ValueError
    else:
        sheet = sheets[0]
    excel_worksheet = excel_workbook[sheet]

    data = []

    i = 2
    a_cell = "A" + str(i)
    b_cell = "B" + str(i)
    c_cell = "C" + str(i)

    while excel_worksheet[a_cell].value is not None and excel_worksheet[b_cell].value is not None:
        if not isinstance(excel_worksheet[a_cell].value, numbers.Number):
            print "ERROR: Value at cell", a_cell, "is not a number"
            raise ValueError
        if not isinstance(excel_worksheet[b_cell].value, numbers.Number):
            print "ERROR: Value at cell", b_cell, "is not a number"
            raise ValueError

        if excel_worksheet[c_cell].value == "Baseline":
            data.append((excel_worksheet[a_cell].value, excel_worksheet[b_cell].value))
        i += 1
        a_cell = "A" + str(i)
        b_cell = "B" + str(i)
        c_cell = "C" + str(i)

    data.sort(key=(lambda x: x[0]))
    time_array = [j[0] for j in data]
    value_array = [j[1] for j in data]

    return np.array(time_array), np.array(value_array)

class Serlfling(object):
    """Serfling sinusoidal regression object."""

    def __init__(self, times, values, cycles):
        self.times = times
        self.values = values

        range = times[-1] - times[0] + 1
        period = range / float(cycles)

        # Regression params: Y(t) = a + b1*t + b2*cos(2*t*pi / 52.1667) + f*sin(2*t*pi / 52.1667) + e
        self.p = period
        self.a = None
        self.b1 = None
        self.b2 = None
        self.b3 = None
        self.e = 0

    def InitializeParameters(self):
        """Initializes the serlfling sinusoidal regression parameters.

        Sets the a, b1, b2, and b3 values of the regression to initial approximations
        suited for optimization.
        """
        lin_reg = linregress(self.times, self.values)
        self.a = lin_reg[1]
        self.b1 = lin_reg[0]
        self.b2 = 3 * np.std(self.values) / (2 ** 0.5)
        self.b3 = 0
        self.e = lin_reg[4]

    def OptimizeParameters(self):
        """Get the fitted alpha, beta 1, and beta 2 parameters.

        Optimizes the current regression parameters according to the input
        time - value pairs.

        Parameters:
            times: (numpy array) the data time values
            values: (numpy array) data values corresponding to each of the time
                entries
        """
        optimize_func = lambda x: (self.a
                                   + self.b1 * self.times
                                   + x[0] * np.sin(self.times * (2 * np.pi / self.p))
                                   + x[1] * np.cos(self.times * (2 * np.pi / self.p))
                                   + self.e
                                   - self.values)
        self.b2, self.b3, = leastsq(optimize_func, [self.b2, self.b3])[0]

    def rmse(self):
        e = (self.a
             + self.b1 * self.times
             + self.b2 * np.sin(self.times * (2 * np.pi / self.p))
             + self.b3 * np.cos(self.times * (2 * np.pi / self.p))
             + self.e
             - self.values)
        return np.sqrt(np.mean((e) ** 2))

    def StdErrorDev(self):
        """Std dev of the error between the data and the regression."""
        e = self.values - (self.a
             + self.b1 * self.times
             + self.b2 * np.sin(self.times * (2 * np.pi / self.p))
             + self.b3 * np.cos(self.times * (2 * np.pi / self.p))
             + self.e
             )

        return np.sqrt(np.mean(e ** 2))

    def PercentileOffset(self, x=95):
        """Calculate offset for x percentile above regreesion.

        Calculates the offset to add to the regression curve for it to be
        above x% of the values in the data set.

        Params:
            x: desired percentile of values below the resulting curve.
        """
        e = self.values - (self.a
             + self.b1 * self.times
             + self.b2 * np.sin(self.times * (2 * np.pi / self.p))
             + self.b3 * np.cos(self.times * (2 * np.pi / self.p))
             + self.e
             )

        e.sort()
        n = len(e)
        return e[int(n * x / 100)]

    def ValuesAbove(self, x = 0):
        n = len(self.values)
        j = 0
        f = self.GetSerflingRegresion()
        for i in xrange(n):
            if f(self.times[i]) + x <= self.values[i]:
                j += 1
        return j



    def GetSerflingRegresion(self):
        """Returns a sinusoidal regression with the current parameters."""
        return lambda t: (self.a
                          + self.b1 * t
                          + self.b2 * np.sin(t * (2 * np.pi / self.p))
                          + self.b3 * np.cos(t * (2 * np.pi / self.p))
                          + self.e)

    def GetPlotValues(self):
        f = self.GetSerflingRegresion()
        data = [f(i) for i in self.times]
        return np.array(data)

    def Print(self):
        print "\n-- Serfling regression data --\n"
        print "period:", self.p
        print "parameters:\n"
        print "a: ", self.a
        print "b1:", self.b1
        print "b2:", self.b2
        print "b3:", self.b3, "\n"
        print "rmse:", self.rmse()
        print "linear regression standard error", self.e, "\n"





# Main ----------

if len(sys.argv[:]) < 3:
    print "ERROR: Not enough arguments."
    print "Make sure to provide the excel file name and the number of cycles in the data"
    raise ValueError

file_name = sys.argv[1]
num_cycles = sys.argv[2]

# We first estimate the regression parameters
print "\nEstimating serfling regression parameters for data in file", file_name, "and", num_cycles, "cycles\n"

time_array, data_array = GetData(file_name)
serfling = Serlfling(time_array, data_array, num_cycles)
n = len(data_array)

serfling.InitializeParameters()
serfling.OptimizeParameters()
serfling.Print()

# Now we estimate useful curve offsets
print ("Estimating percentile offsets and standard deviation of "
       "values above the curve.\n")

percentile95_off = serfling.PercentileOffset(95)
print "- To get 95% of the values below the curve:\n"

print "95% offset:", percentile95_off
print "Values above the curve:", serfling.ValuesAbove(percentile95_off), "of", n
print float(serfling.ValuesAbove(percentile95_off)) / n, "%\n"

percentile975_off = serfling.PercentileOffset(97.5)
print "- To get 97.5% of the values below the curve:\n"

print "97.5% offset:", percentile975_off
print "Values above the curve:", serfling.ValuesAbove(percentile975_off), "of", n
print float(serfling.ValuesAbove(percentile975_off)) / n, "%\n"

std = serfling.StdErrorDev()
p90 = 1.645

print "- Standard deviation of the difference between the curve and the data\n"
print "Std dev:", std, "\n"
print "1.645 std deviation offset:", std * p90
print "Values above the curve:", serfling.ValuesAbove(std * p90), "of", n
print float(serfling.ValuesAbove(std * p90)) / n, "%\n"


p975 = 2


print "2 std deviation offset:", std * p975
print "Values above the curve:", serfling.ValuesAbove(std * p975), "of", n
print float(serfling.ValuesAbove(std * p975)) / n, "%\n"
