from astropy import units as u

HEALPIX_LEVELS = {
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


def coords2hpix(ra, dec, level):
    """
    Input:
     - radec : array shape = (N,2)
     - level : healpix level representing the size of each element
    """
    from healpy.pixelfunc import ang2pix

    nside = 2**level

#     pix_to_visit = []
#     for i in range(len(ra)):
#         ipix = ang2pix(nside, ra[i], dec[i], nest=True, lonlat=True)
#         pix_to_visit.append(ipix)

    def a2p(radec):
        return ang2pix(nside, radec[0], radec[1], nest=True, lonlat=True)

#     pix_to_visit = list(map(a2p, zip(ra,dec)))
    pix_to_visit = map(a2p, zip(ra, dec))

    return pix_to_visit


def hpix2coords(hpixs, level):
    """
    Transform healpix elements to coordinates

    Input:
     - hpixs : list of healpix elements
     - level : healpix level to consider
    """
    from healpy.pixelfunc import pix2ang

    nside = 2**level

#     coords_to_visit = []
#     for ipix in hpixs:
#         c = pixelfunc.pix2ang(nside, ipix, nest=True, lonlat=True)
#         coords_to_visit.append(c)
    def p2a(ipix):
        return pix2ang(nside, ipix, nest=True, lonlat=True)

#     coords_to_visit = list(map(p2a, hpixs))
    coords_to_visit = map(p2a, hpixs)

    return coords_to_visit


def coords_binning_hpix(ra, dec, level, unique=False):
    """
    Return the healpix binning of 'ra,dec' at some 'level'

    Input:
     - radec : array shape = (N,2)
     - level : healpix level representing the size of each element
     - unique: if True, remove (binned) duplicates
    """

    pix_to_visit = coords2hpix(ra, dec, level)

    if unique:
        pix_to_visit = set(pix_to_visit)

    return hpix2coords(pix_to_visit, level)


def coords_binning(ra, dec, rad, unique=False, truncate=False):
    """
    Return list of binned ra,dec in closest healpix size to 'rad'

    Input:
     - ra       : list of (degree) RA coordinates
     - dec      : list of (degree) Dec coordinates
     - rad      : radius/size value to bin;
                  the actual value will be defined by one of Healpix levels;
                  the closest level may be above or below, see 'truncate'
     - unique   : if True, remove duplicated binned coordinates;
                  if False, returned list matches input ra,dec
     - truncate : if True, use the larger-closest healpix element size;
                  if False, use the smaller-closest healpix element
    """
    assert len(ra) == len(dec)
    assert rad > 0

    from moc import utils
    level = utils.size_to_level(rad, truncate=truncate)

    coords_to_visit = coords_binning_hpix(ra, dec, level, unique=unique)

    return coords_to_visit,level

    ## We don't need the neighbours, but if that was the case:
    #
    #nipixs = pixelfunc.get_all_neighbours(nside, c[0], c[1], nest=True, lonlat=True)
    #cns = []
    #for nipix in nipixs:
    #    cns.append(pixelfunc.pix2ang(nside, nipix, nest=True, lonlat=True))
    #coords_to_visit.extend(cns)
