# The apparent position of the Sun at the Vernal Equinox has shifted
# from Aries to Pisces. This script identifies the IAU boundary segment
# that marks this crossing. The coordinates are from the official 1930
# IAU definitions for the B1875.0 epoch.

# The boundary is a line of constant Right Ascension.
ra_h = 23
ra_m = 59
ra_s = 19

# The boundary segment is defined by two declination points.
# Point 1: The southern endpoint.
dec1_deg = -1
dec1_min = 14

# Point 2: The northern endpoint.
dec2_deg = 5
dec2_min = 33

# Convert the declination from degrees and arcminutes to decimal degrees.
# For a negative declination like -1° 14', the value is -(1 + 14/60).
dec1_decimal = dec1_deg - (dec1_min / 60.0)
dec2_decimal = dec2_deg + (dec2_min / 60.0)

# The problem requires the numerically lowest point (by declination) to be first.
# Point 1 has a negative declination, so it comes first.

# Format the Right Ascension string part.
ra_str = f"{ra_h} {ra_m} {ra_s}"

# Format each point's string according to "XX YY ZZ, AA.BB".
# The f-string format ':.2f' rounds the declination to two decimal places.
point1_str = f"{ra_str}, {dec1_decimal:.2f}"
point2_str = f"{ra_str}, {dec2_decimal:.2f}"

# Combine the two points into the final answer string, separated by a semicolon.
final_answer = f"{point1_str}; {point2_str}"

print(final_answer)