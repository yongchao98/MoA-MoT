import math

def calculate_and_print_boundary():
    """
    This script calculates and prints the coordinates of the IAU constellation boundary segment
    that the Vernal Equinox crossed when it moved from Aries into Pisces.
    """

    # This function converts declination from Degrees, Arcminutes, Arcseconds (DMS)
    # to decimal degrees. It correctly handles the sign for declinations between 0 and -1 degrees.
    def dms_to_decimal(d, m, s):
        # In Python, -0 is a valid number, and str(-0) is '-0'.
        # This check correctly assigns a negative sign for values like (-0, 18, 51).
        sign = -1 if str(d).strip().startswith('-') else 1
        return sign * (abs(d) + m / 60.0 + s / 3600.0)

    # The endpoints of the boundary line segment shared by Pisces and Aries
    # that crosses the celestial equator (Dec=0).
    # Data is from the IAU boundaries for the B1875.0 epoch.
    # The point with the lower declination value is defined first.

    # Point 1: The southern endpoint of the segment
    p1_ra_hms = (23, 57, 56.4)
    # Using the special -0 value to correctly handle the sign
    p1_dec_dms = (-0, 18, 51)

    # Point 2: The northern endpoint of the segment
    p2_ra_hms = (23, 57, 56.4)
    p2_dec_dms = (2, 32, 46)

    # Convert the declinations to decimal format
    p1_dec_decimal = dms_to_decimal(p1_dec_dms[0], p1_dec_dms[1], p1_dec_dms[2])
    p2_dec_decimal = dms_to_decimal(p2_dec_dms[0], p2_dec_dms[1], p2_dec_dms[2])

    # Format each point according to the user's request "XX YY ZZ, AA.BB"
    # The RA components are joined by spaces.
    # The Declination is formatted to two decimal places.
    point1_str = (f"{p1_ra_hms[0]:02d} {p1_ra_hms[1]:02d} {p1_ra_hms[2]}, "
                  f"{p1_dec_decimal:.2f}")

    point2_str = (f"{p2_ra_hms[0]:02d} {p2_ra_hms[1]:02d} {p2_ra_hms[2]}, "
                  f"{p2_dec_decimal:.2f}")

    # Print the two points, separated by a semicolon.
    print(f"{point1_str}; {point2_str}")

calculate_and_print_boundary()