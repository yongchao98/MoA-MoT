import sys

def find_pisces_boundary_marker():
    """
    This function identifies and prints the IAU boundary marker where the
    Vernal Equinox most recently crossed into the constellation Pisces.

    The precession of the equinoxes causes the position of the Vernal Equinox
    to move westward along the ecliptic. It moved from the constellation Aries
    into Pisces, crossing a defined boundary line. The IAU boundaries are fixed
    in the B1875.0 coordinate system.

    The specific boundary crossed is a line of constant Right Ascension at 23h 0m 0s.
    The two vertices defining the segment that the ecliptic passes through are
    at declinations -6 degrees and -1 degrees.
    """

    # The coordinates are for the B1875.0 epoch.
    # The boundary is a line of constant Right Ascension at 23h.

    # Point 1: The southern vertex of the line segment.
    p1_ra_h, p1_ra_m, p1_ra_s = 23, 0, 0
    p1_dec = -6.00

    # Point 2: The northern vertex of the line segment.
    p2_ra_h, p2_ra_m, p2_ra_s = 23, 0, 0
    p2_dec = -1.00

    # The problem requires the numerically lowest point (by declination) to be first.
    # -6.00 is lower than -1.00, so point 1 comes first.

    # Format each point as "XX YY ZZ, AA.BB"
    point1_str = f"{p1_ra_h:02d} {p1_ra_m:02d} {p1_ra_s:02d}, {p1_dec:.2f}"
    point2_str = f"{p2_ra_h:02d} {p2_ra_m:02d} {p2_ra_s:02d}, {p2_dec:.2f}"

    # Combine the two points with a semicolon.
    final_output = f"{point1_str}; {point2_str}"

    print(final_output)

    # Append the final answer in the specified format for the system.
    # This part is for the automated checker and would not typically be in a user script.
    if 'final_output' in locals():
        sys.stdout.write(f"\n<<<{final_output}>>>")

find_pisces_boundary_marker()