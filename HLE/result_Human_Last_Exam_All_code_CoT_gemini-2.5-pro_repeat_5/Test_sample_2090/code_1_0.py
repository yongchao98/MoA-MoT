def find_constellation_boundary_marker():
    """
    This function identifies and formats the specific IAU constellation boundary
    marker where the Vernal Equinox point crossed from Aries into Pisces.

    The boundary is a line segment at a constant Right Ascension. The data is
    based on the official IAU constellation boundaries defined by E. Delporte in 1930.
    """

    # The boundary segment is at a constant Right Ascension of 00h 40m 00s.
    ra_h = 0
    ra_m = 40
    ra_s = 0

    # The segment is defined by two declination points.
    # Point 1: +02 degrees 30 minutes 00 seconds
    # Point 2: +06 degrees 30 minutes 00 seconds
    dec1_deg = 2.0
    dec1_min = 30.0

    dec2_deg = 6.0
    dec2_min = 30.0

    # Convert declinations to decimal degrees for comparison and formatting.
    p1_dec_decimal = dec1_deg + dec1_min / 60.0
    p2_dec_decimal = dec2_deg + dec2_min / 60.0

    # Ensure the point with the numerically lowest declination comes first.
    # In this case, +2.50 is lower than +6.50.
    
    # Format the two points into the required string format "XX YY ZZ, AA.BB".
    # The format specifier {:02d} ensures leading zeros for single-digit numbers.
    # The format specifier {:.2f} ensures two decimal places for the declination.
    point1_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {p1_dec_decimal:.2f}"
    point2_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {p2_dec_decimal:.2f}"

    # Combine the two points with a semicolon.
    final_answer = f"{point1_str}; {point2_str}"
    
    print(final_answer)

find_constellation_boundary_marker()