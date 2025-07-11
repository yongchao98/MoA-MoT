def find_pisces_boundary_marker():
    """
    This script identifies and formats the IAU constellation boundary segment
    that the Vernal Equinox crossed when moving into Pisces.
    """
    # The boundary is between Cetus and Pisces. The vertices are defined
    # in the B1875.0 epoch coordinate system.

    # Point 1: The southern vertex of the segment.
    ra1_h, ra1_m, ra1_s = 0, 0, 0
    dec1_deg, dec1_min = -7, 30

    # Point 2: The northern vertex of the segment.
    ra2_h, ra2_m, ra2_s = 0, 0, 0
    dec2_deg, dec2_min = 6, 0

    # Convert declination to decimal degrees.
    # For the southern point, both degrees and minutes are negative.
    dec1_decimal = dec1_deg - (dec1_min / 60.0)
    # For the northern point.
    dec2_decimal = dec2_deg + (dec2_min / 60.0)

    # Format the points as "XX YY ZZ, AA.BB"
    # The numerically lowest point (by declination, since RAs are equal) comes first.
    point1_formatted = f"{ra1_h:02d} {ra1_m:02d} {ra1_s:02d}, {dec1_decimal:.2f}"
    point2_formatted = f"{ra2_h:02d} {ra2_m:02d} {ra2_s:02d}, {dec2_decimal:.2f}"

    # Combine the two points with a semicolon as requested.
    final_answer = f"{point1_formatted}; {point2_formatted}"

    print(final_answer)

find_pisces_boundary_marker()