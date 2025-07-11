import math

def format_boundary_points():
    """
    This function defines the IAU boundary points, formats them according to the specified
    rules, and prints the final result.
    The boundary is where the Vernal Equinox most recently crossed into Pisces.
    This corresponds to the IAU B1875.0 epoch boundary between Pisces and Cetus.
    """

    # The boundary is a line of constant Right Ascension.
    # The two points defining the segment that crosses the celestial equator are:
    # Point A: RA = 00h 06m 40.0s, Dec = +06deg 30min
    # Point B: RA = 00h 06m 40.0s, Dec = -02deg 30min

    # Store RA components (same for both points)
    ra_h = 0
    ra_m = 6
    ra_s = 40

    # Store Declination in decimal degrees
    # +06deg 30min = 6.5 degrees
    # -02deg 30min = -2.5 degrees
    dec_high = 6.5
    dec_low = -2.5

    # The prompt asks for the numerically lowest point first.
    # Since RAs are identical, we order by declination.
    point1_dec = dec_low
    point2_dec = dec_high

    # Format the points into the "XX YY ZZ, AA.BB" format.
    # Point 1 (lower declination)
    p1_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {point1_dec:.2f}"

    # Point 2 (higher declination)
    p2_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {point2_dec:.2f}"

    # Combine the two formatted strings with a semicolon.
    final_output = f"{p1_str}; {p2_str}"

    print(final_output)

format_boundary_points()