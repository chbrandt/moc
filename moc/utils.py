import itertools
import numpy as np
from astropy import units
from .core import HEALPIX_LEVELS


def size_to_level(size, truncate=False):
    '''
    Return Healpix level matching to 'size'

    The original idea was to define the size of Healpix diamonds
    based on the positional error size. As a default, the safe
    size of a diamond has been chosen to be twice the size.

    Input:
     - size : float or astropy.Quantity
        Size to consider for defining Healpix level.
        If not an astropy.Quantity, 'size' is assumed to be in "degree"
     - truncate : bool
        If True, return the closest, but greater-than-size level;
        Otherwise, return the closest, but smaller-than-size level

    Output:
     - Healpix level : int
        Returned value correspond to the closest HEALpix size to 'size'
    '''
    try:
        size.unit
    except AttributeError as e:
        size = size * units.degree

    levels = list(HEALPIX_LEVELS.keys())
    levels.sort()

    ko = None
    for i in levels:
        ko = i
        hpx_size = HEALPIX_LEVELS[i]
        if hpx_size < size:
            break
    ko = max(0, ko - int(truncate))
    return ko


def naive_region2pointings(ramin, ramax, decmin, decmax, radius):

    bbox = dict(ramin=ramin * units.degree,
                ramax=ramax * units.degree,
                decmin=decmin * units.degree,
                decmax=decmax * units.degree)

    try:
        rad = radius.to('arcmin')
    except:
        rad = radius * units.arcmin

    # get the closest smaller size of healpix elements
    level = size_to_level(rad, truncate=True)
    #
    dsize = HEALPIX_LEVELS[level]

    # create a grid of (fake) coordinates, to then create a MOC table from it
    step_size = dsize.to('deg').value/2**0.5

    ra_vec = np.arange(bbox['ramin'].value, bbox['ramax'].value, step_size)
    dec_vec = np.arange(bbox['decmin'].value, bbox['decmax'].value, step_size)

    grid = np.asarray(list(itertools.product(ra_vec, dec_vec))).T

    return grid
