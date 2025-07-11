import math

def solve():
    """
    This script determines and formats the IAU constellation boundary
    that the Vernal Equinox crossed when moving from Aries to Pisces.
    """

    # The IAU boundary between Aries and Pisces is a line of constant
    # Right Ascension at 0h 0m 0s (for the B1875.0 epoch).
    # The segment crossing the ecliptic is defined by two vertices.

    # Coordinates for the first point (southern endpoint)
    p1_ra_h, p1_ra_m, p1_ra_s = 0, 0, 0
    p1_dec_d, p1_dec_m, p1_dec_s = -6, 5, 0

    # Coordinates for the second point (northern endpoint)
    p2_ra_h, p2_ra_m, p2_ra_s = 0, 0, 0
    p2_dec_d, p2_dec_m, p2_dec_s = 5, 55, 0

    # Convert declination to decimal degrees.
    # For negative declinations, both degrees and minutes contribute negatively.
    p1_dec_decimal = p1_dec_d - (p1_dec_m / 60.0) - (p1_dec_s / 3600.0)
    p2_dec_decimal = p2_dec_d + (p2_dec_m / 60.0) + (p2_dec_s / 3600.0)

    # The problem asks for the numerically lowest point to be first.
    # Since the RAs are identical, we order by declination.
    # Point 1 has a declination of -6.08, Point 2 has a declination of 5.92.
    # So, Point 1 comes first.

    # Format each point into the required string: "XX YY ZZ, AA.BB"
    point1_str = f"{p1_ra_h:02d} {p1_ra_m:02d} {p1_ra_s:02d}, {p1_dec_decimal:.2f}"
    point2_str = f"{p2_ra_h:02d} {p2_ra_m:02d} {p2_ra_s:02d}, {p2_dec_decimal:.2f}"

    # Combine the two points with a semicolon as a separator.
    final_boundary_line = f"{point1_str}; {point2_str}"

    print(final_boundary_line)

solve()
<<<00 00 00, -6.08; 00 00 00, 5.92>>>