import math

def convert_dec_to_decimal(degrees, minutes, seconds):
    """
    Converts a Declination from degrees, minutes, seconds to decimal degrees.
    Handles negative declinations correctly.
    """
    sign = -1 if degrees < 0 else 1
    decimal_degrees = abs(degrees) + minutes / 60.0 + seconds / 3600.0
    return sign * decimal_degrees

# The B1875.0 coordinates for the vertices of the boundary segment
# between Aries and Pisces that crosses the celestial equator.
# This segment is a line of constant Right Ascension.
# Source: Delporte, E. (1930), Délimitation scientifique des constellations.

# Point 1 (southern point)
ra1_h, ra1_m, ra1_s = 23, 59, 14
dec1_d, dec1_m, dec1_s = -10, 4, 31

# Point 2 (northern point)
ra2_h, ra2_m, ra2_s = 23, 59, 14
dec2_d, dec2_m, dec2_s = 14, 55, 29

# Convert declinations to decimal format
dec1_decimal = convert_dec_to_decimal(dec1_d, dec1_m, dec1_s)
dec2_decimal = convert_dec_to_decimal(dec2_d, dec2_m, dec2_s)

# Format the points as per the required specification "XX YY ZZ, AA.BB"
# The numerically lowest point (Point 1 with negative declination) comes first.
# Ensure integer seconds for RA are displayed as integers.
point1_str = f"{ra1_h} {ra1_m} {ra1_s}, {dec1_decimal:.2f}"
point2_str = f"{ra2_h} {ra2_m} {ra2_s}, {dec2_decimal:.2f}"

# Combine the two points into the final answer string
final_answer = f"{point1_str}; {point2_str}"

# Print the final result
print("The boundary marker is defined by the line segment connecting the two points:")
print(final_answer)
