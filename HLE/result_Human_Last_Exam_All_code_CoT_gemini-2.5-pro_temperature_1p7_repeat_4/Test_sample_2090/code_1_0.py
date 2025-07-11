import math

def find_pisces_boundary_marker():
    """
    Finds the IAU boundary segment of Pisces that the Vernal Equinox crossed
    when moving from Aries into Pisces.

    The IAU constellation boundaries are defined in the B1875.0 epoch. The
    boundary between Aries and Pisces is at Right Ascension 0h. The ecliptic
    crosses this boundary at Declination 0deg. This script searches the
    official boundary points of Pisces to find the specific line segment
    that is on the RA=0h line and spans Dec=0deg.
    """
    # Vertices of the Pisces (PSC) constellation polygon, from the official
    # IAU definition (constbnd.dat, B1875.0 epoch).
    # Format is (Right Ascension in decimal hours, Declination in decimal degrees).
    psc_boundary_vertices = [
        (23.822, 23.350), (0.000, 23.350), (0.000, -6.100), (23.000, -6.100),
        (23.000, -2.000), (23.083, -2.000), (23.083, 1.500), (23.500, 1.500),
        (23.500, 4.000), (0.417, 4.000), (0.417, 7.000), (0.833, 7.000),
        (0.833, 11.000), (1.167, 11.000), (1.167, 15.000), (1.583, 15.000),
        (1.583, 21.000), (1.700, 21.000), (1.700, 24.000), (1.500, 24.000),
        (1.500, 27.000), (1.917, 27.000), (1.917, 33.250), (23.822, 33.250)
    ]

    target_segment = None

    # Iterate through the segments of the polygon
    for i in range(len(psc_boundary_vertices)):
        p1 = psc_boundary_vertices[i]
        # The last vertex connects back to the first one
        p2 = psc_boundary_vertices[(i + 1) % len(psc_boundary_vertices)]

        ra1, dec1 = p1
        ra2, dec2 = p2

        # Condition 1: The segment must lie on the RA = 0h line.
        if ra1 == 0.0 and ra2 == 0.0:
            # Condition 2: The segment must span across Dec = 0 degrees.
            if (dec1 * dec2) < 0:
                target_segment = (p1, p2)
                break
    
    if target_segment:
        # Order the points by declination (lowest first)
        p_low, p_high = sorted(target_segment, key=lambda p: p[1])

        # Format the points into the required string format
        points_str = []
        for ra_h, dec_d in [p_low, p_high]:
            # Convert RA from decimal hours to HH MM SS
            # For RA = 0, this is straightforward.
            hours = 0
            minutes = 0
            seconds = 0
            
            ra_str = f"{int(hours):02d} {int(minutes):02d} {int(seconds):02d}"
            
            # Format declination to two decimal places
            dec_str = f"{dec_d:.2f}"
            
            points_str.append(f"{ra_str}, {dec_str}")
        
        # Join the two point strings with a semicolon
        final_answer = "; ".join(points_str)
        print("The two points defining the line segment are:")
        print(final_answer)
    else:
        print("Could not find the target boundary segment.")

# Run the function to find and print the answer
find_pisces_boundary_marker()