# This script prints the coordinates of the IAU constellation boundary segment
# where the Vernal Equinox point crossed from Aries into Pisces.

# The boundary segment is a line of constant Right Ascension.
# The coordinates are based on the B1875.0 epoch as defined by the IAU.

# Coordinates for the first point of the line segment
ra_h = 1
ra_m = 44
ra_s = 0
dec_1 = 10.00  # Declination in decimal degrees

# Coordinates for the second point of the line segment
# Right Ascension is the same.
dec_2 = 15.00  # Declination in decimal degrees

# Format the output string for each point according to "XX YY ZZ, AA.BB"
# The point with the lower declination is presented first.
point_1_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {dec_1:.2f}"
point_2_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {dec_2:.2f}"

# Combine the two points into the final answer string, separated by a semicolon.
final_answer = f"{point_1_str}; {point_2_str}"

print(final_answer)