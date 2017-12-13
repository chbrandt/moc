from astropy.table import Table
from astropy import units
from mocpy import MOC

from .utils import size_to_level


def table_to_moc(table, ra, dec, radius, factor=2, outfile=None):
    '''
    Return a MOC table from 'table['ra','dec']'

    Input:
     - table : astropy.Table
        Table containing columns 'ra' and 'dec'
     - ra : string
        Name of 'table' RA column
     - dec : string
        Name of 'table' Dec column
     - radius : float or astropy.Quantity
        Together with 'factor' defines the size of HEALPix diamonds.
        Typically this will be the (mean) positional error.
        If not a astropy.Quantity, 'radius' is considered to be in 'deg'
     - factor : float
        Together with 'radius' defines the size of HEALPix diamonds.
        For instance, a "factor=2" means the region's diameter
     - outfile : string
        Optional name for the output MOC catalog to be written to

    Output:
     - MOC catalog
    '''
    assert isinstance(table, Table)

    try:
        radius.unit
    except AttributeError as e:
        radius = radius * units.degree

    size = factor*radius
    level = size_to_level(size)

    moc = MOC.from_table(table, ra, dec, level)
    if outfile:
        moc.write(outfile)
    return moc
