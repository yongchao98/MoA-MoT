import sys

def find_pisces_boundary_marker():
    """
    This function provides the specific IAU boundary marker where the Vernal Equinox
    crossed into the constellation Pisces.

    The boundary between Aries and Pisces is a line of constant Right Ascension
    defined in the B1875.0 epoch. The two points defining this line segment are
    retrieved from the official IAU constellation boundary definitions.
    """

    # The boundary is defined by two points with the same Right Ascension.
    # The point with the lower declination value is listed first.

    # Point 1: The southern endpoint, at the junction of Pisces, Cetus, and Aries.
    # Right Ascension: 0 hours, 0 minutes, 0 seconds
    # Declination: -9.50 degrees
    ra1_h, ra1_m, ra1_s = 0, 0, 0
    dec1 = -9.50

    # Point 2: The northern endpoint, at the junction of Pisces, Andromeda, and Aries.
    # Right Ascension: 0 hours, 0 minutes, 0 seconds
    # Declination: 6.00 degrees
    ra2_h, ra2_m, ra2_s = 0, 0, 0
    dec2 = 6.00

    # Format the coordinates into the specified "XX YY ZZ, AA.BB" string format.
    # Using format specifiers to ensure leading zeros for alignment.
    point1_str = f"{ra1_h:02d} {ra1_m:02d} {ra1_s:02d}, {dec1:.2f}"
    point2_str = f"{ra2_h:02d} {ra2_m:02d} {ra2_s:02d}, {dec2:.2f}"

    # Combine the two points, separated by a semicolon, as requested.
    final_boundary_marker = f"{point1_str}; {point2_str}"

    print(final_boundary_marker)
    # Use sys.stdout.flush() to ensure the output is printed before the marker.
    sys.stdout.flush() 

if __name__ == "__main__":
    find_pisces_boundary_marker()
