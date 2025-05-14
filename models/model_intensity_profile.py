import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from ALMA_specs import alma_dict


def gaussian_ring_profile(R, r1, sigma1):
    """
    Return a 1D Gaussian ring intensity profile.
    Parameters:
        R : ndarray
            Array of radial positions (same units as r1 and sigma1)
        r1 : float
            Radius of the ring (AU)
        sigma1 : float
            Width (standard deviation) of the ring (AU)
    Returns:
        I(R) : ndarray
            Intensity profile evaluated at each radius R (dimensionless)
    """
    return np.exp(-((R - r1)**2.) / (2. * sigma1**2.))


def radial_to_image(profile_r, r_array, npix=512, pixel_scale=1.0, inclination=0.0, PA=0.0):
    """
    Convert a radial profile into a 2D image, applying inclination and position angle.
    Parameters:
        profile_r : ndarray
            Intensity as a function of radius (same length as r_array)
        r_array : ndarray
            Radial positions corresponding to profile_r (in AU)
        npix : int
            Number of pixels per side of the output image
        pixel_scale : float
            Size of a pixel in same units as r_array (e.g., AU)
        inclination : float
            Inclination angle in degrees (0 = face-on, 90 = edge-on)
        PA : float
            Position angle in degrees (counter-clockwise from North)
    Returns:
        image : 2D ndarray
            Inclined, rotated 2D image of the disk
    """
    # Shift origin to center
    y, x = np.indices((npix, npix))
    center = (npix - 1) / 2.
    x = (x - center) * pixel_scale
    y = (y - center) * pixel_scale

    # Rotate coordinates by -PA (counter-clockwise, from North to East)
    pa_rad = np.radians(PA)
    x_rot = x * np.cos(pa_rad) + y * np.sin(pa_rad)
    y_rot = -x * np.sin(pa_rad) + y * np.cos(pa_rad)

    # Apply inclination (squash x-axis)
    inc_rad = np.radians(inclination)
    x_rot /= np.cos(inc_rad)

    # Compute projected radius
    r_grid = np.sqrt(x_rot**2 + y_rot**2)

    # Interpolate and reshape
    image = np.interp(r_grid.ravel(), r_array, profile_r)
    return image.reshape((npix, npix))



def save_image_to_fits(image, pixel_scale_au, output_path, config, npix):
    """
    Save a 2D image to a FITS file with full WCS header based on config.
    Parameters:
        image : 2D ndarray
            The image to be saved
        pixel_scale_au : float
            Pixel size in AU
        output_path : str
            Path to output FITS file
        config : dict
            Dictionary with simulation parameters (RA, Dec, distance, freq...)
        npix : int
            Image size (number of pixels per side)
    """
    distance_pc = config["distance"]
    pixel_scale_arcsec = (pixel_scale_au / distance_pc) 
    px_size_deg = pixel_scale_arcsec / 3600

    band = config["ALMA_band"]
    specs = alma_dict[band]
    central_freq = config.get("central_freq", specs["central_freq"])

    # Convert RA/Dec to degrees
    coord = SkyCoord(config["RA"], config["Dec"], unit=(u.hourangle, u.deg))
    ra_deg = coord.ra.deg
    dec_deg = coord.dec.deg

    # Create FITS header
    hdr = fits.Header()
    hdr['BTYPE'] = 'Intensity'
    hdr['OBJECT'] = 'Disk model'
    hdr['BUNIT'] = 'Jy/pixel'

    # Spatial axes
    hdr['CTYPE1'] = 'RA---SIN'
    hdr['CRVAL1'] = ra_deg
    hdr['CDELT1'] = -px_size_deg
    hdr['CRPIX1'] = npix / 2
    hdr['CUNIT1'] = 'deg'

    hdr['CTYPE2'] = 'DEC--SIN'
    hdr['CRVAL2'] = dec_deg
    hdr['CDELT2'] = px_size_deg
    hdr['CRPIX2'] = npix / 2
    hdr['CUNIT2'] = 'deg'

    # Frequency axis
    hdr['CTYPE3'] = 'FREQ'
    hdr['CRVAL3'] = central_freq
    hdr['CDELT3'] = -1e5  # placeholder value for frequency resolution
    hdr['CRPIX3'] = 1.0
    hdr['CUNIT3'] = 'Hz'

    # Stokes axis
    hdr['CTYPE4'] = 'STOKES'
    hdr['CRVAL4'] = 1.0
    hdr['CDELT4'] = 1.0
    hdr['CRPIX4'] = 1.0
    hdr['CUNIT4'] = ''

    # Additional metadata
    hdr['RESTFRQ'] = central_freq
    hdr['SPECSYS'] = 'LSRK'

    # Write FITS
    hdu = fits.PrimaryHDU(data=image, header=hdr)
    hdul = fits.HDUList([hdu])
    hdul.writeto(output_path, overwrite=True)
    print("Saved FITS image to:", output_path)
