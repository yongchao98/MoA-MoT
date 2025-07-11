import math

def solve_boundary_crossing():
    """
    Finds the IAU constellation boundary segment that the Vernal Equinox
    most recently crossed into Pisces.
    """

    # Vertices for the eastern boundary of Pisces bordering Aries, derived from
    # the IAU definition (VizieR catalog VI/49A, J2000.0).
    # The format is (RA_H, RA_M, Dec_D, Dec_M).
    boundary_points_data = [
        (1, 57.8, 3, 18),
        (1, 57.8, 7, 3),
        (2, 6.9, 7, 3),
        (2, 6.9, 10, 39),
        (2, 8.4, 10, 39),
        (2, 8.4, 15, 47),
    ]

    # Convert raw data into a list of dictionaries with decimal degrees/hours
    vertices = []
    for p in boundary_points_data:
        ra_h, ra_m, dec_d, dec_m = p
        ra_decimal_h = ra_h + ra_m / 60.0
        dec_decimal_d = dec_d + dec_m / 60.0
        vertices.append({'ra_h': ra_decimal_h, 'dec_d': dec_decimal_d})

    # Obliquity of the ecliptic for J2000.0 in degrees
    epsilon_deg = 23.4392911
    epsilon_rad = math.radians(epsilon_deg)
    tan_epsilon = math.tan(epsilon_rad)

    intersecting_segment = None

    # Iterate through segments defined by consecutive vertices
    for i in range(len(vertices) - 1):
        p1 = vertices[i]
        p2 = vertices[i+1]

        # Check vertical segments (constant RA)
        if p1['ra_h'] == p2['ra_h']:
            ra_const_deg = p1['ra_h'] * 15
            ra_const_rad = math.radians(ra_const_deg)
            
            # Find the declination of the ecliptic at this RA
            dec_ecliptic_rad = math.atan(tan_epsilon * math.sin(ra_const_rad))
            dec_ecliptic_deg = math.degrees(dec_ecliptic_rad)

            # Check if the intersection point is on the segment
            dec_min = min(p1['dec_d'], p2['dec_d'])
            dec_max = max(p1['dec_d'], p2['dec_d'])
            
            if dec_min <= dec_ecliptic_deg <= dec_max:
                intersecting_segment = (p1, p2)
                # Since we process in increasing RA, this will be the one.
                break 

    if not intersecting_segment:
        print("Error: Could not find the intersecting boundary segment.")
        return

    # Format the coordinates of the two points of the segment
    def format_point(point):
        # RA to HH MM SS
        ra_h_total = point['ra_h']
        h = int(ra_h_total)
        m_rem = (ra_h_total - h) * 60
        m = int(m_rem)
        s_rem = (m_rem - m) * 60
        s = int(round(s_rem))
        ra_str = f"{h:02d} {m:02d} {s:02d}"

        # Dec to AA.BB
        dec_d_total = point['dec_d']
        dec_str = f"{dec_d_total:.2f}"
        
        return f"{ra_str}, {dec_str}"

    p1, p2 = intersecting_segment

    # Order points by declination, lowest first
    if p1['dec_d'] > p2['dec_d']:
        p1, p2 = p2, p1
    
    point1_str = format_point(p1)
    point2_str = format_point(p2)

    # Per the instructions, outputting each number in the final equation
    p1_ra_formatted = point1_str.split(', ')[0]
    p1_dec_formatted = point1_str.split(', ')[1]
    p2_ra_formatted = point2_str.split(', ')[0]
    p2_dec_formatted = point2_str.split(', ')[1]

    print("The boundary marker is a line segment defined by two points.")
    print(f"Point 1 (RA, Dec): {p1_ra_formatted}, {p1_dec_formatted}")
    print(f"Point 2 (RA, Dec): {p2_ra_formatted}, {p2_dec_formatted}")

    # Combine into the final answer string
    final_answer = f"{point1_str}; {point2_str}"
    print("\nFormatted Answer:")
    print(final_answer)


solve_boundary_crossing()
<<<02 08 24, 10.65; 02 08 24, 15.78>>>