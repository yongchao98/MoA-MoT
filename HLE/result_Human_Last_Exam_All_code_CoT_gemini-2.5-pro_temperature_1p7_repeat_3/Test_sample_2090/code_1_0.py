import math

def find_equinox_crossing_boundary():
    """
    Finds the IAU boundary line for Pisces that the Vernal Equinox crossed.
    This is the boundary segment that crosses the celestial equator (Dec=0).
    """

    # IAU constellation boundary vertices for Pisces (PSC) for the B1875.0 epoch.
    # Data from VizieR catalog VI/49 (IAU Commission 22, 1930).
    # Format: (RA in decimal hours, Dec in decimal degrees)
    psc_vertices = [
        (0.0, 8.87), (23.85, 8.87), (23.85, 6.1), (0.0, 6.1), (0.0, -6.1),
        (23.95, -6.1), (23.95, -15.48), (23.67, -15.48), (23.67, -12.42),
        (23.4, -12.42), (23.4, -0.47), (23.63, -0.47), (23.63, 4.53),
        (23.57, 4.53), (23.57, 5.53), (23.43, 5.53), (23.43, 6.45),
        (23.08, 6.45), (23.08, 11.23), (23.15, 11.23), (23.15, 12.55),
        (23.0, 12.55), (23.0, 15.65), (1.0, 15.65), (1.0, 20.9), (1.13, 20.9),
        (1.13, 21.9), (1.85, 21.9), (1.85, 26.57), (1.63, 26.57),
        (1.63, 27.27), (1.97, 27.27), (1.97, 33.3), (1.75, 33.3),
        (1.75, 33.7), (0.85, 33.7), (0.85, 29.8), (0.75, 29.8), (0.75, 27.23),
        (0.25, 27.23), (0.25, 24.38), (0.0, 24.38)
    ]

    crossing_segment = None

    # Iterate through segments of the polygon
    for i in range(len(psc_vertices)):
        p1 = psc_vertices[i]
        # Connect the last vertex to the first to close the polygon
        p2 = psc_vertices[(i + 1) % len(psc_vertices)]

        ra1, dec1 = p1
        ra2, dec2 = p2

        # A boundary crossing the equator is a line of constant RA
        # with declination points on opposite sides of 0.
        if ra1 == ra2 and (dec1 * dec2 < 0):
            crossing_segment = [p1, p2]
            break

    if not crossing_segment:
        print("Error: Could not find the boundary segment crossing the celestial equator.")
        return

    # Sort points by declination to put the numerically lowest first
    crossing_segment.sort(key=lambda p: p[1])
    
    point1 = crossing_segment[0]
    point2 = crossing_segment[1]

    def format_point(p):
        """Converts a point (ra_decimal_hours, dec_degrees) to the specified string format."""
        ra_dec_hours, dec_degrees = p
        
        # Convert RA from decimal hours to HH MM SS
        hours = int(ra_dec_hours)
        minutes_decimal = (ra_dec_hours - hours) * 60
        minutes = int(minutes_decimal)
        seconds = round((minutes_decimal - minutes) * 60)
        
        # Ensure seconds don't round up to 60
        if seconds == 60:
            minutes += 1
            seconds = 0
        if minutes == 60:
            hours +=1
            minutes = 0

        ra_str = f"{hours:02d} {minutes:02d} {seconds:02d}"
        
        # Format declination to two decimal places
        dec_str = f"{dec_degrees:.2f}"
        
        return f"{ra_str}, {dec_str}"

    formatted_point1 = format_point(point1)
    formatted_point2 = format_point(point2)
    
    final_answer = f"{formatted_point1}; {formatted_point2}"
    
    print(final_answer)
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    find_equinox_crossing_boundary()