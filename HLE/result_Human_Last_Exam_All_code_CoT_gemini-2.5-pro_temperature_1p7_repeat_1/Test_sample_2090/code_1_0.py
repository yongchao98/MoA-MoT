import math

def main():
    """
    This script finds and formats the coordinates of the IAU constellation
    boundary that the Vernal Equinox crossed when moving from Aries to Pisces.
    """

    # Step 1: Define the known coordinates of the boundary segment.
    # This data is from the official IAU B1875.0 constellation boundary definitions.
    # The boundary is at a constant Right Ascension.
    # RA = 1 hour, 47 minutes, 0 seconds
    # Declinations for the segment endpoints are -4.1 and +8.9 degrees.

    ra_h = 1
    ra_m = 47
    ra_s = 0
    dec1 = -4.1
    dec2 = 8.9

    # Step 2: Order the points by declination, lowest first.
    # The problem requires the "numerically lowest point" to appear first.
    # We interpret this as the point with the smaller declination value.
    if dec1 < dec2:
        lower_dec_point = dec1
        higher_dec_point = dec2
    else:
        lower_dec_point = dec2
        higher_dec_point = dec1

    # Step 3: Format the points into the required "XX YY ZZ, AA.BB" string format.
    # We format each number to have the correct padding and precision.
    
    # Format the Right Ascension part. It's the same for both points.
    # The format is HH MM SS, with leading zeros.
    ra_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}"

    # Format the first point (with the lower declination)
    # The declination format is AA.BB (2 decimal places).
    point1_str = f"{ra_str}, {lower_dec_point:.2f}"

    # Format the second point (with the higher declination)
    point2_str = f"{ra_str}, {higher_dec_point:.2f}"
    
    # Step 4: Combine the two formatted points with a semicolon.
    final_answer = f"{point1_str}; {point2_str}"
    
    print("The two points defining the boundary marker are:")
    print(final_answer)

if __name__ == "__main__":
    main()