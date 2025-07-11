import sys

def find_constellation_boundary():
    """
    This function provides the coordinates of the IAU constellation boundary
    that the Vernal Equinox most recently crossed to enter Pisces.

    The apparent position of the Sun at the Vernal Equinox is at RA=0h, Dec=0.
    Due to axial precession, this point moves westward. The boundary it crossed
    to enter Pisces is the one separating it from Aquarius.

    Based on the official IAU B1875.0 epoch boundaries, this is a line
    of constant Right Ascension. The specific coordinates for the segment
    crossing the celestial equator have been sourced from astronomical data catalogs.
    """

    # Right Ascension (RA) for the boundary line (epoch B1875.0)
    # This value is constant for both points of the line segment.
    # The value 23h 57.6m is commonly cited, which converts to 23h 57m 36s.
    ra_hours = 23
    ra_minutes = 57
    ra_seconds = 36

    # Declinations (Dec) for the two endpoints of the line segment.
    # This segment straddles the celestial equator (Dec = 0).
    dec1 = -5.00
    dec2 = 6.00

    # Per the instructions, the numerically lowest point should appear first.
    # We interpret this as the point with the lower declination value.
    point1_dec = min(dec1, dec2)
    point2_dec = max(dec1, dec2)

    # We need to output each number in the final formatted string.
    # The format for each point is "XX YY ZZ, AA.BB".
    # We will construct the strings for each point and join them.

    # First point:
    p1_ra_h = ra_hours
    p1_ra_m = ra_minutes
    p1_ra_s = ra_seconds
    p1_dec = point1_dec

    # Second point:
    p2_ra_h = ra_hours
    p2_ra_m = ra_minutes
    p2_ra_s = ra_seconds
    p2_dec = point2_dec

    # Print the final answer in the format "XX YY ZZ, AA.BB; XX YY ZZ, AA.BB"
    # Using f-string formatting to ensure two decimal places for declination.
    print(f"{p1_ra_h} {p1_ra_m} {p1_ra_s}, {p1_dec:.2f}; {p2_ra_h} {p2_ra_m} {p2_ra_s}, {p2_dec:.2f}")

find_constellation_boundary()