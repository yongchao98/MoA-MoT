import math

def solve_boundary_question():
    """
    This function identifies and prints the specific IAU constellation boundary
    that the Vernal Equinox crossed when moving into Pisces.
    """

    # The IAU B1875.0 constellation boundaries for Pisces (Psc) define the
    # border with Aries and Cetus along the Right Ascension = 0h line.
    # The ecliptic (and the celestial equator at the equinox) crosses this
    # line between the following two vertices.

    # Point 1: RA = 0h 0m 0s, Dec = -2.8 degrees
    # Point 2: RA = 0h 0m 0s, Dec = +6.5 degrees
    
    boundary_points = [
        {'ra_h': 0, 'ra_m': 0, 'ra_s': 0.0, 'dec': -2.8},
        {'ra_h': 0, 'ra_m': 0, 'ra_s': 0.0, 'dec': 6.5}
    ]

    # The problem asks for the numerically lowest point to appear first.
    # We will sort by declination.
    boundary_points.sort(key=lambda p: p['dec'])

    # Format the points into the required string format: "XX YY ZZ, AA.BB"
    output_parts = []
    for point in boundary_points:
        ra_h = point['ra_h']
        ra_m = point['ra_m']
        # Use math.trunc to handle the integer part of seconds for formatting
        ra_s = math.trunc(point['ra_s'])
        dec = point['dec']
        
        # Format: RA as 0-padded integers, Dec as a float with 2 decimal places
        formatted_point = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {dec:.2f}"
        output_parts.append(formatted_point)

    # Join the formatted points with a semicolon
    final_answer = "; ".join(output_parts)
    
    print("The boundary marker is defined by the line segment connecting the following two points:")
    print(final_answer)

solve_boundary_question()