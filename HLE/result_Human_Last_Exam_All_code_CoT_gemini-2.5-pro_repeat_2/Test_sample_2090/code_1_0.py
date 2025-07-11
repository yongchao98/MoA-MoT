import math

def solve():
    """
    Finds and formats the IAU constellation boundary segment that the Vernal Equinox
    crossed when moving into Pisces.
    """
    # According to the IAU definitions (B1875.0 epoch), the boundary between
    # Pisces and Aries/Cetus is along the line of Right Ascension = 0h.
    # We list the vertices along this line.
    # Format: (RA_hours, Declination_degrees)
    # The declination values are converted from degrees and arcminutes.
    # +5d 45m = 5 + 45/60 = 5.75
    # -10d 30m = -10 - 30/60 = -10.5
    boundary_vertices = [
        (0.0, 20.75),  # Corner with Andromeda
        (0.0, 5.75),   # Corner with Aries
        (0.0, -10.5),  # Corner with Cetus
        (0.0, -25.0)   # Southern corner
    ]

    # The Vernal Equinox is at Declination = 0. We need to find the boundary
    # segment that crosses the celestial equator (Dec=0).
    target_segment = None
    for i in range(len(boundary_vertices) - 1):
        p1 = boundary_vertices[i]
        p2 = boundary_vertices[i+1]
        
        # A segment crosses Dec=0 if one endpoint has a positive declination
        # and the other has a negative declination.
        if p1[1] * p2[1] < 0:
            target_segment = (p1, p2)
            break
            
    if not target_segment:
        print("Error: Could not find the target boundary segment.")
        return

    # The problem asks for the "numerically lowest point" to appear first.
    # Since RA is the same for both points, we sort by declination.
    p_low, p_high = sorted(target_segment, key=lambda p: p[1])

    # This function formats a point (RA_h, Dec_deg) into the required string.
    def format_point(point):
        ra_h, dec_d = point
        
        # Convert RA from decimal hours to H, M, S.
        # Since input RA is 0.0, the result is 0, 0, 0.
        frac_h, h = math.modf(ra_h)
        frac_m, m = math.modf(frac_h * 60)
        s = frac_m * 60
        
        ra_hours = int(h)
        ra_minutes = int(m)
        ra_seconds = int(round(s))
        
        # Format the numbers into the final strings as per the requirement.
        # RA: "XX YY ZZ", Dec: "AA.BB"
        ra_str = f"{ra_hours:02d} {ra_minutes:02d} {ra_seconds:02d}"
        dec_str = f"{dec_d:.2f}"
        
        return f"{ra_str}, {dec_str}"

    # Format the two points of the segment.
    point1_str = format_point(p_low)
    point2_str = format_point(p_high)
    
    # The final answer is the two points separated by a semicolon.
    final_output = f"{point1_str}; {point2_str}"
    
    print(final_output)

solve()
<<<00 00 00, -10.50; 00 00 00, 5.75>>>