"""
uitl.ncfunc : A set of useful functions when dealing with netcdf4 data.

This module provides functions to help deal with netCDF data. 

Functions list:
    * get_nc_file(fname, rw)
    * copy_att(nc_source, nc_target)
    * copy_atts_bad_fill(nc_source, nc_target, missing_value)
    * copy_atts_add_fill(nc_source, nc_target, missing_value)
"""

from netCDF4 import Dataset
import numpy as np


def get_nc_file(fname, rw):
    """Get a netcdf file.
    """
    nc_file = Dataset(fname, rw, format="NETCDF4")
    return nc_file


def copy_atts(nc_source, nc_target):
    """Copy netCDF attributes.

    This function copies the attributes from one netCDF element to another.
    
    Parameters
    ----------
    nc_source :
        Source netCDF element
    nc_target :
        Target netCDF element

    Examples
    --------
    Copy the attributes from one variable to another.

    >>> old_var = nc_source.variables['old']
    >>> new_var = nc_target.createVariable('new', 'f4', ('y','x',) )
    >>> new_var[:,:] = old_var[:,:]
    >>> copy_atts( old_var,new_var )
    """

    # get a list of global attribute names from the incoming file
    atts = nc_source.ncattrs()

    # place those attributes in the outgoing file
    for ii in range(len(atts)):
        nc_target.setncattr(atts[ii], nc_source.getncattr(atts[ii]))


def copy_atts_bad_fill(nc_source, nc_target, missing_value):
    """Copy all netCDF attributes except _FillValue.  

    This function copies all the attributes from one netCDF element to another,
    but ignores the _FillValue attribute and sets MissingValue. 
    
    Parameters
    ----------
    nc_source :
        Source netCDF element
    nc_target :
        Target netCDF element
    missing_value :
        Value to set as indicator of missing values.

    Examples
    --------
    Copy the attributes from one variable to another.

    >>> old_var = nc_source.variables['old']
    >>> new_var = nc_target.createVariable('new', 'f4', ('y','x',) )
    >>> new_var[:,:] = old_var[:,:]
    >>> copy_atts_bad_fill( old_var,new_var, -9999. )
    """

    # get a list of global attribute names from the incoming file
    atts = nc_source.ncattrs()

    # place those attributes in the outgoing file
    for ii in range(len(atts)):
        if atts[ii] != "_FillValue":
            nc_target.setncattr(atts[ii], nc_source.getncattr(atts[ii]))
        else:
            nc_target.setncattr("missing_value", np.float32(missing_value))


def copy_atts_add_fill(nc_source, nc_target, missing_value):
    """Copy all netCDF attributes and add a missing value attribute.  

    This function copies all the attributes from one netCDF element to another 
    and adds a missing value attribute to nc_target.
    
    Parameters
    ----------
    nc_source :
        Source netCDF element
    nc_target :
        Target netCDF element
    missing_value :
        Value to set as indicator of missing values.

    Examples
    --------
    Copy the attributes from one variable to another.

    >>> old_var = nc_source.variables['old']
    >>> new_var = nc_target.createVariable('new', 'f4', ('y','x',) )
    >>> new_var[:,:] = old_var[:,:]
    >>> copy_atts_bad_fill( old_var,new_var, -9999. )
    """

    # get a list of global attribute names from the incoming file
    atts = nc_source.ncattrs()

    # place those attributes in the outgoing file
    for ii in range(len(atts)):
        nc_target.setncattr(atts[ii], nc_source.getncattr(atts[ii]))

    nc_target.setncattr("missing_value", np.float32(missing_value))
