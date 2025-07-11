import sys
import subprocess

# Ensure astropy is installed
try:
    from astropy.coordinates import SkyCoord
    from astropy import units as u
except ImportError:
    print("Astropy library not found. Please install it using: pip install astropy")
    sys.exit(1)

def find_pisces_boundary():
    """
    This function calculates the modern J2000.0 coordinates for the IAU boundary
    segment where the ecliptic enters Pisces from Aquarius.
    """

    # The IAU constellation boundaries are defined in the B1875.0 coordinate system.
    # The relevant boundary segment separating Aquarius and Pisces that is crossed
    # by the ecliptic is a line of constant Right Ascension.
    # Its endpoints in B1875.0 are:
    p1_ra_b1875 = '23h09m03.1s'
    p1_dec_b1875 = '-06d45m00s'

    p2_ra_b1875 = '23h09m03.1s'
    p2_dec_b1875 = '-09d15m00s'

    # Create Astropy SkyCoord objects for these points using the FK4 frame for B1875.0
    try:
        point1_b1875 = SkyCoord(ra=p1_ra_b1875, dec=p1_dec_b1875, frame='fk4', equinox='B1875.0')
        point2_b1875 = SkyCoord(ra=p2_ra_b1875, dec=p2_dec_b1875, frame='fk4', equinox='B1875.0')
    except Exception as e:
        print(f"An error occurred creating SkyCoord objects: {e}")
        return

    # Transform the coordinates to the modern J2000.0 system (ICRS)
    point1_j2000 = point1_b1875.transform_to('icrs')
    point2_j2000 = point2_b1875.transform_to('icrs')

    # Order the points by Right Ascension, which is the standard convention.
    if point1_j2000.ra.value < point2_j2000.ra.value:
        first_point = point1_j2000
        second_point = point2_j2000
    else:
        first_point = point2_j2000
        second_point = point1_j2000

    # Format the coordinates as specified: RA as HH MM SS, Dec as DD.DD
    # First point
    fp_ra_hms = first_point.ra.hms
    fp_ra_str = f"{int(fp_ra_hms.h):02d} {int(fp_ra_hms.m):02d} {int(round(fp_ra_hms.s)):02d}"
    fp_dec_str = f"{first_point.dec.degree:.2f}"
    first_point_formatted = f"{fp_ra_str}, {fp_dec_str}"

    # Second point
    sp_ra_hms = second_point.ra.hms
    sp_ra_str = f"{int(sp_ra_hms.h):02d} {int(sp_ra_hms.m):02d} {int(round(sp_ra_hms.s)):02d}"
    sp_dec_str = f"{second_point.dec.degree:.2f}"
    second_point_formatted = f"{sp_ra_str}, {sp_dec_str}"

    # Print the final result
    print(f"{first_point_formatted}; {second_point_formatted}")

if __name__ == '__main__':
    find_pisces_boundary()
<<<23 32 17, -6.51; 23 32 27, -8.97>>>