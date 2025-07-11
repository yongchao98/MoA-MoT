import math

# The IAU constellation boundary between Pisces and its neighbors (Aries/Cetus)
# that crosses the celestial equator is a line of constant Right Ascension.
# These values are for the B1875.0 epoch.

# Right Ascension for the boundary line
ra_hours = 23
ra_minutes = 49
ra_seconds = 19

# The two points defining the line segment that crosses the celestial equator (Dec = 0)
# Point 1 has a negative declination.
# Point 2 has a positive declination.
point1_dec = -6.26
point2_dec = 5.53

# The problem asks for the numerically lowest point first.
# Since the Right Ascension is the same, this refers to the declination.
# Point with declination -6.26 comes before the point with declination 5.53.

# Format the points into the required string "XX YY ZZ, AA.BB"
point1_str = f"{ra_hours:02d} {ra_minutes:02d} {ra_seconds:02d}, {point1_dec:.2f}"
point2_str = f"{ra_hours:02d} {ra_minutes:02d} {ra_seconds:02d}, {point2_dec:.2f}"

# Combine the two points with a semicolon as a separator.
final_answer = f"{point1_str}; {point2_str}"

# Print the final formatted string.
# This prints all the numbers that constitute the final answer as requested.
print(final_answer)