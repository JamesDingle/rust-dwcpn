import xarray as xr
import numpy

import sys


filename=sys.argv[1]
var = sys.argv[2]

ds = xr.open_dataset(filename)

var_data = ds[var]

pp_masked = numpy.ma.masked_where(numpy.isnan(var_data), var_data)

pp_no_masked_vals = numpy.asarray(pp_masked[~pp_masked.mask])

print(f"Mean: {pp_no_masked_vals.mean()}")
print(f"Median: {numpy.median(pp_no_masked_vals)}")
print(f"1st percentile: {numpy.percentile(pp_no_masked_vals, 1)}")
print(f"25th percentile: {numpy.percentile(pp_no_masked_vals, 25)}")
print(f"75th percentile: {numpy.percentile(pp_no_masked_vals, 75)}")
print(f"99th percentile: {numpy.percentile(pp_no_masked_vals, 99)}")
# print(f"Mean: {pp.mean()}")