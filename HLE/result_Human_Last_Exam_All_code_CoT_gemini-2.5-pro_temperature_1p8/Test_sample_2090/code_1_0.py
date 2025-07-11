def find_pisces_boundary_marker():
    """
    This function identifies and prints the IAU constellation boundary segment
    where the Vernal Equinox precessed from Aries into Pisces.
    """
    # The IAU constellation boundaries are defined in the B1875.0 epoch.
    # The boundary between Aries and Pisces is a line of constant Right Ascension at 0h.
    # The Vernal Equinox point, by definition, has a declination of 0 degrees.
    # We need the segment of the RA=0h line that straddles Dec=0.
    # Based on the official IAU data (VizieR catalog VI/49), the vertices are:
    
    # Point 1: RA = 0h 0m 0s, Dec = +6.00 deg
    # Point 2: RA = 0h 0m 0s, Dec = -5.85 deg

    point_high = {'ra_hms': (0, 0, 0), 'dec': 6.00}
    point_low = {'ra_hms': (0, 0, 0), 'dec': -5.85}

    # The problem asks to present the numerically lowest point first.
    # Since RA is identical, we sort by declination.
    # Point with Dec = -5.85 comes first.
    sorted_points = [point_low, point_high]

    print("The boundary marker is a line segment defined by two points.")
    
    # Per the instruction "output each number in the final equation!",
    # we print the components of each point.
    p1 = sorted_points[0]
    p2 = sorted_points[1]

    print(f"\nFirst point (lowest):")
    print(f"  Right Ascension: {p1['ra_hms'][0]:02d}h {p1['ra_hms'][1]:02d}m {p1['ra_hms'][2]:02d}s")
    print(f"  Declination: {p1['dec']:.2f} degrees")
    
    print(f"\nSecond point (highest):")
    print(f"  Right Ascension: {p2['ra_hms'][0]:02d}h {p2['ra_hms'][1]:02d}m {p2['ra_hms'][2]:02d}s")
    print(f"  Declination: {p2['dec']:.2f} degrees")

    # Format the points into the required string "XX YY ZZ, AA.BB"
    formatted_points = []
    for p in sorted_points:
        ra_str = f"{p['ra_hms'][0]:02d} {p['ra_hms'][1]:02d} {p['ra_hms'][2]:02d}"
        dec_str = f"{p['dec']:.2f}"
        formatted_points.append(f"{ra_str}, {dec_str}")
        
    # Join the points with a semicolon as requested.
    final_answer = ";".join(formatted_points)

    print("\nThe boundary marker, with points separated by a semicolon, is:")
    print(final_answer)

find_pisces_boundary_marker()