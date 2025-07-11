def find_constellation_boundary():
    """
    This function determines and prints the IAU constellation boundary segment
    that the Vernal Equinox most recently crossed into Pisces.
    """

    # The boundary between Aries and Pisces is at Right Ascension 0h.
    # The Vernal Equinox is at (RA=0h, Dec=0), so we need the boundary
    # segment that crosses Dec=0.
    # According to the IAU B1875.0 boundary definitions, the two vertices
    # on the RA=0h line that bracket the celestial equator are:
    # Point 1: RA = 0h 0m 0s, Dec = -12 deg
    # Point 2: RA = 0h 0m 0s, Dec = +6 deg

    # Define the coordinates for the first point (lower declination)
    p1_ra_h = 0
    p1_ra_m = 0
    p1_ra_s = 0
    p1_dec = -12.0

    # Define the coordinates for the second point (higher declination)
    p2_ra_h = 0
    p2_ra_m = 0
    p2_ra_s = 0
    p2_dec = 6.0

    # Format each point according to "XX YY ZZ, AA.BB"
    # The f-string formatting ensures leading zeros for single-digit numbers
    # and two decimal places for the declination.
    point1_str = f"{p1_ra_h:02d} {p1_ra_m:02d} {p1_ra_s:02d}, {p1_dec:.2f}"
    point2_str = f"{p2_ra_h:02d} {p2_ra_m:02d} {p2_ra_s:02d}, {p2_dec:.2f}"

    # Combine the two points with a semicolon, with the numerically
    # lowest (by declination) point first.
    final_boundary_string = f"{point1_str}; {point2_str}"

    print(final_boundary_string)

find_constellation_boundary()