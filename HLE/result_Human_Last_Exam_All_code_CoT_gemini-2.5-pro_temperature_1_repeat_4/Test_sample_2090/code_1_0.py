import math

def solve_constellation_boundary():
    """
    This script finds the specific IAU constellation boundary segment that the
    vernal equinox point crossed into Pisces.
    """
    
    # Relevant constellation boundary data for Pisces (PSC) from the official
    # IAU B1875.0 definitions (CDS catalog VI/42).
    # Format: RA (HH MM.m), Dec (DD MM), Constellation
    raw_data = """
    23 45.8 -06 32 PSC
    00 00.0 -06 32 PSC
    00 00.0 +02 05 PSC
    00 26.8 +02 05 PSC
    """

    # Parse the raw data into a list of points
    points = []
    for line in raw_data.strip().split('\n'):
        parts = line.strip().split()
        if not parts:
            continue
        
        ra_h, ra_m_float = int(parts[0]), float(parts[1])
        dec_d_str, dec_m_str = parts[2], parts[3]
        
        # Store original components for precise final formatting
        point_data = {
            "ra_comp": (ra_h, ra_m_float),
            "dec_comp": (int(dec_d_str), int(dec_m_str))
        }
        points.append(point_data)

    # Find the boundary segment at RA=0h that crosses the celestial equator
    target_segment = None
    for i in range(len(points) - 1):
        p1 = points[i]
        p2 = points[i+1]

        # Check for a line of constant RA at 0h
        is_ra_zero_line = (p1["ra_comp"][0] == 0 and p1["ra_comp"][1] == 0.0 and
                           p2["ra_comp"][0] == 0 and p2["ra_comp"][1] == 0.0)

        if is_ra_zero_line:
            # Calculate declination values to check for equator crossing
            dec1_val = p1["dec_comp"][0] + math.copysign(p1["dec_comp"][1] / 60.0, p1["dec_comp"][0])
            dec2_val = p2["dec_comp"][0] + math.copysign(p2["dec_comp"][1] / 60.0, p2["dec_comp"][0])

            # Check if segment crosses Dec=0
            if dec1_val * dec2_val < 0:
                # Ensure the point with lower declination is first
                if dec1_val < dec2_val:
                    target_segment = (p1, p2)
                else:
                    target_segment = (p2, p1)
                break

    if not target_segment:
        print("Could not find the specified boundary segment.")
        return

    # Format the two points for the final output
    formatted_points = []
    for point in target_segment:
        # Format RA: XX YY ZZ
        ra_h, ra_m_float = point["ra_comp"]
        ra_m = int(ra_m_float)
        ra_s = int(round((ra_m_float - ra_m) * 60.0))
        ra_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}"

        # Format Dec: AA.BB
        dec_d, dec_m = point["dec_comp"]
        # Calculate decimal degrees from original components for accuracy
        dec_val = dec_d + math.copysign(dec_m / 60.0, dec_d if dec_d != 0 else 1.0)
        dec_str = f"{dec_val:+.2f}"
        
        formatted_points.append(f"{ra_str}, {dec_str}")
    
    # Print the final result separated by a semicolon
    print("; ".join(formatted_points))

solve_constellation_boundary()