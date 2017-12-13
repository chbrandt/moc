from astropy import units as u

healpix_levels = {
    0  : 58.63 * u.deg,
    1  : 29.32 * u.deg,
    2  : 14.66 * u.deg,
    3  : 7.329 * u.deg,
    4  : 3.665 * u.deg,
    5  : 1.832 * u.deg,
    6  : 54.97 * u.arcmin,
    7  : 27.48 * u.arcmin,
    8  : 13.74 * u.arcmin,
    9  : 6.871 * u.arcmin,
    10 : 3.435 * u.arcmin,
    11 : 1.718 * u.arcmin,
    12 : 51.53 * u.arcsec,
    13 : 25.77 * u.arcsec,
    14 : 12.88 * u.arcsec,
    15 : 6.442 * u.arcsec,
    16 : 3.221 * u.arcsec,
    17 : 1.61 * u.arcsec,
    18 : 805.2 * u.milliarcsecond,
    19 : 402.6 * u.milliarcsecond,
    20 : 201.3 * u.milliarcsecond,
    21 : 100.6 * u.milliarcsecond,
    22 : 50.32 * u.milliarcsecond,
    23 : 25.16 * u.milliarcsecond,
    24 : 12.58 * u.milliarcsecond,
    25 : 6.291 * u.milliarcsecond,
    26 : 3.145 * u.milliarcsecond,
    27 : 1.573 * u.milliarcsecond,
    28 : 786.3 * u.microarcsecond,
    29 : 393.2 * u.microarcsecond
}

def radius_to_level(pos_error,factor=2):
    '''
    Return Healpix level matching to 'radius'

    The original idea was to define the size of Healpix diamonds
    based on the positional error radius. As a default, the safe
    size of a diamond has been chosen to be twice the error radius.

    Input:
     - pos_error : float
        Radius to consider for defining Healpix level
     - factor : float
        Multiplicative value over 'pos_error'

    Output:
     - Healpix level : int
        Returned value correspond to the closest
        size greater-or-equal to "factor * pos_error"
    '''
    assert pos_error.unit
    ko = 0
    levels = healpix_levels.keys()
    levels.sort()
    for lvl in levels:
        rad = healpix_levels[lvl]
        if rad < factor * pos_error:
            break
        ko = lvl
    return ko

def radec_to_moc(table, ra, dec, pos_error, output_filename):
    from astropy.table import Table
    assert isinstance(table,Table)

    level = radius_to_level(pos_error)

    from mocpy import MOC
    moc = MOC.from_table(table,ra,dec,level)
    moc.write(output_filename)
