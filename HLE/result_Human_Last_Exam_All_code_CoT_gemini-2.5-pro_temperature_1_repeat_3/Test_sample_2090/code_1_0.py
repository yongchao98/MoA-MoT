import sys

def find_equinox_crossing_boundary():
    """
    This script finds the IAU constellation boundary segment that the Vernal Equinox
    crossed when it moved into Pisces.
    
    The IAU constellation boundaries are defined in the B1875.0 epoch. The boundary
    between Aries and Pisces is located at Right Ascension 0h. The Vernal Equinox
    is located on the celestial equator (Declination 0 deg). Therefore, we need
    to find the line segment at RA=0h that crosses Dec=0 deg.
    """

    # These are the official IAU B1875.0 boundary vertices for Pisces along the
    # 0h Right Ascension line, which it shares with Aries and Cetus.
    # Data format: (RA_hours, RA_minutes, RA_seconds, Declination_degrees)
    pisces_ra0_vertices = [
        (0, 0, 0, 32.50),
        (0, 0, 0, 27.25),
        (0, 0, 0, 2.50),
        (0, 0, 0, -6.50),
        (0, 0, 0, -12.75)
    ]

    crossing_segment = None
    # Iterate through pairs of adjacent vertices to find the segment crossing Dec=0
    for i in range(len(pisces_ra0_vertices) - 1):
        p1 = pisces_ra0_vertices[i]
        p2 = pisces_ra0_vertices[i+1]
        
        # A segment crosses the equator if one point's declination is positive
        # and the other's is negative.
        if p1[3] * p2[3] < 0:
            # Check if it's on the RA=0 line
            if p1[0] == 0 and p1[1] == 0 and p1[2] == 0:
                 crossing_segment = (p1, p2)
                 break

    if not crossing_segment:
        print("Error: Could not find the crossing boundary segment.", file=sys.stderr)
        return

    # The problem asks to list the numerically lowest point first.
    # Since RA is the same, this refers to the declination.
    sorted_segment = sorted(crossing_segment, key=lambda p: p[3])
    
    p_low = sorted_segment[0]
    p_high = sorted_segment[1]

    # Extract each number for the final equation/string.
    # Point 1 (lowest declination)
    rah1, ram1, ras1, dec1 = p_low
    # Point 2 (highest declination)
    rah2, ram2, ras2, dec2 = p_high

    print("The boundary marker is the line segment connecting two points.")
    print("The numbers defining the first point (lowest declination) are:")
    print(f"RA: {rah1}h {ram1}m {ras1}s, Declination: {dec1:.2f} degrees")
    print("\nThe numbers defining the second point are:")
    print(f"RA: {rah2}h {ram2}m {ras2}s, Declination: {dec2:.2f} degrees")

    # Format the points into the required string format: XX YY ZZ, AA.BB
    point1_str = f"{rah1:02d} {ram1:02d} {ras1:02d}, {dec1:.2f}"
    point2_str = f"{rah2:02d} {ram2:02d} {ras2:02d}, {dec2:.2f}"

    # Combine into the final answer, separated by a semicolon
    final_answer = f"{point1_str}; {point2_str}"
    
    print("\nFinal formatted answer string:")
    print(final_answer)

find_equinox_crossing_boundary()