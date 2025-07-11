import sys

def find_boundary_marker():
    """
    This script identifies and prints the coordinates of the IAU constellation boundary line
    that the vernal equinox crossed when moving from Aries to Pisces.

    The boundary is a line of constant Right Ascension (0h 0m 0s) and is defined by two
    vertices that bracket the celestial equator.
    """

    # Coordinates of the southern vertex of the boundary line segment
    ra_h_south = 0
    ra_m_south = 0
    ra_s_south = 0
    dec_south = -9.00

    # Coordinates of the northern vertex of the boundary line segment
    ra_h_north = 0
    ra_m_north = 0
    ra_s_north = 0
    dec_north = 6.00

    # The user has requested to output each number in the final equation.
    # We will print the components of the two points that define the line segment.
    # The lowest point (by declination) should appear first.
    
    # Printing components for the first point (southern, lower declination)
    print("Point 1 (South) Right Ascension (HH MM SS):", f"{ra_h_south:02d}", f"{ra_m_south:02d}", f"{ra_s_south:02d}", file=sys.stdout)
    print("Point 1 (South) Declination (DD.DD):", f"{dec_south:.2f}", file=sys.stdout)
    
    # Printing components for the second point (northern, higher declination)
    print("Point 2 (North) Right Ascension (HH MM SS):", f"{ra_h_north:02d}", f"{ra_m_north:02d}", f"{ra_s_north:02d}", file=sys.stdout)
    print("Point 2 (North) Declination (DD.DD):", f"{dec_north:.2f}", file=sys.stdout)
    
    # Format the points into the final string as requested.
    point1_str = f"{ra_h_south:02d} {ra_m_south:02d} {ra_s_south:02d}, {dec_south:.2f}"
    point2_str = f"{ra_h_north:02d} {ra_m_north:02d} {ra_s_north:02d}, {dec_north:.2f}"

    # Combine the two points with a semicolon for the final answer format.
    final_answer = f"{point1_str}; {point2_str}"
    
    print("\nFinal Answer Format (Point1; Point2):")
    print(final_answer)

find_boundary_marker()

# The final answer is then presented in the required format.
# The string represents the two endpoints of the line segment.
final_formatted_answer = "00 00 00, -9.00; 00 00 00, 6.00"
print(f'<<<{final_formatted_answer}>>>')