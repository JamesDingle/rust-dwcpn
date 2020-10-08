import xarray as xr
import numpy
filename="pp_processing_20100501_30.000_50.000_-140.000_-120.000.nc"

ds = xr.open_dataset(filename)

pp = ds['pp']

pp_masked = numpy.ma.masked_where(numpy.isnan(pp), pp)

pp_no_masked_vals = numpy.asarray(pp_masked[~pp_masked.mask])

print(f"Mean: {pp_no_masked_vals.mean()}")
print(f"Median: {numpy.median(pp_no_masked_vals)}")
print(f"1st percentile: {numpy.percentile(pp_no_masked_vals, 1)}")
print(f"25th percentile: {numpy.percentile(pp_no_masked_vals, 25)}")
print(f"75th percentile: {numpy.percentile(pp_no_masked_vals, 75)}")
print(f"99th percentile: {numpy.percentile(pp_no_masked_vals, 99)}")
# print(f"Mean: {pp.mean()}")