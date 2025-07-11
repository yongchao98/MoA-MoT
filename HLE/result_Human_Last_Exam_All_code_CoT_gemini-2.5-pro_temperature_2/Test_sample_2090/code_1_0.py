import math

def format_declination(dec_degrees):
    """Formats the declination to two decimal places, handling positive signs."""
    return f"{dec_degrees:.2f}"

def main():
    """
    This program prints the coordinates of the IAU constellation boundary
    segment that separates Pisces and Aries. This is the line the
    Vernal Equinox point crossed most recently to enter Pisces.
    The coordinates are for the B1875.0 epoch, as defined by the IAU.
    """
    # The boundary between Aries and Pisces is a line of constant Right Ascension.
    # Data is from the original IAU (Delporte, 1930) definitions.
    ra_hours = 23
    ra_minutes = 22
    ra_seconds = 30

    # The line segment is defined by two declination points.
    # Point 1 (Southern vertex, on the boundary with Cetus)
    dec1_degrees = -6.75  # -6 degrees and 45 arcminutes

    # Point 2 (Northern vertex, on the boundary with Andromeda)
    dec2_degrees = 14.0   # +14 degrees and 0 arcminutes

    # The problem asks for the numerically lowest point first.
    # Since the RA is the same, we compare by declination.
    # Point 1 has a lower declination than Point 2.

    # Format the points as XX YY ZZ, AA.BB
    point1_str = f"{ra_hours:02d} {ra_minutes:02d} {ra_seconds:02d}, {format_declination(dec1_degrees)}"
    point2_str = f"{ra_hours:02d} {ra_minutes:02d} {ra_seconds:02d}, {format_declination(dec2_degrees)}"

    # Print the final result separated by a semicolon.
    final_output = f"{point1_str}; {point2_str}"
    print(final_output)

if __name__ == "__main__":
    main()