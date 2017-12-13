from .core import HEALPIX_LEVELS


def size_to_level(size, truncate=False):
    '''
    Return Healpix level matching to 'size'

    The original idea was to define the size of Healpix diamonds
    based on the positional error size. As a default, the safe
    size of a diamond has been chosen to be twice the size.

    Input:
     - size : float
        size to consider for defining Healpix level
     - truncate : bool
        If True, return the closest, but greater-than-size level;
        Otherwise, return the closest, but smaller-than-size level

    Output:
     - Healpix level : int
        Returned value correspond to the closest HEALpix size to 'size'
    '''
    assert size.unit
    ko = None
    levels = HEALPIX_LEVELS.keys()
    levels.sort()
    for i in levels:
        ko = i
        hpx_size = HEALPIX_LEVELS[i]
        if hpx_size < size:
            break
    ko = max(0, ko - int(truncate))
    return ko
