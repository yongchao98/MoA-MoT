def solve_constellation_boundary():
    """
    This script formats and prints the IAU constellation boundary
    where the Vernal Equinox point entered Pisces from Aries.
    """

    # The boundary is a line segment at RA = 0h (B1875.0 epoch).
    # The endpoints are defined by their declination.
    # Point 1: RA 0h 0m 0s, Dec +5d 0m 0s
    # Point 2: RA 0h 0m 0s, Dec -6d 30m 0s

    points_data = [
        {'ra': (0, 0, 0), 'dec': (5, 0, 0)},
        {'ra': (0, 0, 0), 'dec': (-6, 30, 0)}
    ]

    # Helper function to convert declination to decimal degrees
    def to_decimal_degrees(dec_tuple):
        d, m, s = dec_tuple
        sign = -1 if d < 0 else 1
        return d + sign * (m / 60.0) + sign * (s / 3600.0)

    # Process and format each point
    formatted_points = []
    for p in points_data:
        ra_h, ra_m, ra_s = p['ra']
        
        # Get declination as a single decimal value for sorting and printing
        dec_decimal = to_decimal_degrees(p['dec'])
        
        # Format the point according to the specified format "XX YY ZZ, AA.BB"
        point_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {dec_decimal:.2f}"
        
        formatted_points.append({
            'dec': dec_decimal,
            'str': point_str
        })

    # Sort the points by declination, lowest first
    formatted_points.sort(key=lambda x: x['dec'])
    
    # Join the formatted strings with a semicolon
    final_answer = "; ".join([p['str'] for p in formatted_points])
    
    print(final_answer)

solve_constellation_boundary()