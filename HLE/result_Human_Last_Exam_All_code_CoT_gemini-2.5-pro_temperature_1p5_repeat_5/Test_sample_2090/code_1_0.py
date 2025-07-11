# The apparent position of the Sun at the Vernal Equinox has shifted over
# the years due to precession. This script identifies the official IAU
# boundary line that this position crossed when it moved from the
# constellation Aries into Pisces.
#
# This boundary, for the B1875.0 epoch, is a line of constant Right
# Ascension. The specific segment that crosses the celestial equator
# (and therefore the ecliptic at the equinox) is what we need.

# Coordinates of the boundary line segment between Aries and Pisces
# that crosses the celestial equator.
# Source: IAU 1930 definitions (Delporte, B1875.0 epoch).

# Right Ascension (constant for both points)
ra_hours = 23
ra_minutes = 26
ra_seconds = 19 # Rounded from 19.1s for the required format

# Declination for the two endpoints of the segment
dec_point1_degrees = -5.00
dec_point2_degrees = 5.00

# The problem asks for the numerically lowest point first.
# In this case, that refers to the declination.

# Format the first point (lower declination)
# Format: XX YY ZZ, AA.BB
# We use f-string formatting to ensure leading zeros where needed.
# e.g., {dec_point1_degrees:05.2f} means format as a float with 2 decimal
# places, a total width of 5 characters, padded with a leading zero if needed.
# This handles the sign as well, so -5.00 becomes -05.00.
point1_str = f"{ra_hours:02d} {ra_minutes:02d} {ra_seconds:02d}, {dec_point1_degrees:05.2f}"

# Format the second point (higher declination)
point2_str = f"{ra_hours:02d} {ra_minutes:02d} {ra_seconds:02d}, {dec_point2_degrees:05.2f}"

# Combine the two points with a semicolon as the separator
final_boundary_str = f"{point1_str}; {point2_str}"

# Print the final formatted string
print(final_boundary_str)