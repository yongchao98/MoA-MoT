import math

def format_ra(hours, minutes, seconds):
    """Formats Right Ascension into XX YY ZZ string."""
    return f"{int(hours):02d} {int(minutes):02d} {int(seconds):02d}"

def format_dec(degrees):
    """Formats Declination into AA.BB string."""
    return f"{degrees:.2f}"

def solve_pisces_boundary():
    """
    This function determines and prints the IAU boundary marker
    where the Vernal Equinox crossed into Pisces.
    """
    # The boundary is a line of constant Right Ascension at 0h 0m 0s.
    ra_h, ra_m, ra_s = 0, 0, 0

    # The endpoints of this boundary segment that crosses the celestial equator
    # are at Declination -5 degrees and +20 degrees.
    dec1_deg = -5.0
    dec2_deg = 20.0

    # Ensure the points are ordered by the numerically lowest declination first.
    points_dec = sorted([dec1_deg, dec2_deg])

    # Format the two points.
    # Point 1 (southernmost)
    ra_str_p1 = format_ra(ra_h, ra_m, ra_s)
    dec_str_p1 = format_dec(points_dec[0])
    point1_str = f"{ra_str_p1}, {dec_str_p1}"

    # Point 2 (northernmost)
    ra_str_p2 = format_ra(ra_h, ra_m, ra_s)
    dec_str_p2 = format_dec(points_dec[1])
    point2_str = f"{ra_str_p2}, {dec_str_p2}"

    # Combine the two points into the final answer string.
    final_answer = f"{point1_str}; {point2_str}"

    print(final_answer)

solve_pisces_boundary()
<<<00 00 00, -5.00; 00 00 00, 20.00>>>