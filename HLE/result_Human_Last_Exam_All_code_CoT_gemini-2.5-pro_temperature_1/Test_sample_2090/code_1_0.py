import math

# Plan:
# The IAU boundary between the constellations of Aries and Pisces is a line of constant
# Right Ascension at 0 hours in the B1875.0 coordinate system. The Vernal Equinox,
# which is the point where the ecliptic crosses the celestial equator (Declination = 0),
# crossed this line. We need to identify the two IAU-defined vertices on this line
# that bracket the celestial equator.

# According to the official IAU boundary definitions, these two points are:
# Point 1 (Southern vertex): RA = 0h 0m 0s, Dec = -7 degrees 30 arcminutes
# Point 2 (Northern vertex): RA = 0h 0m 0s, Dec = +6 degrees 0 arcminutes

# The script will format these two points as requested and print the result.

# Define the coordinates for the two vertices of the boundary segment.
# The "numerically lowest point" is interpreted as the one with the lower declination.

# Vertex 1 (numerically lowest declination)
ra1_h, ra1_m, ra1_s = 0, 0, 0
dec1_deg, dec1_min = -7, 30

# Vertex 2 (higher declination)
ra2_h, ra2_m, ra2_s = 0, 0, 0
dec2_deg, dec2_min = 6, 0

# Convert declinations to decimal degrees.
# For dec1, the value is negative, so we subtract the minutes part.
dec1_decimal = dec1_deg - (dec1_min / 60.0)
dec2_decimal = dec2_deg + (dec2_min / 60.0)

# Format each number into the required "XX YY ZZ, AA.BB" format.
# The instruction "output each number in the final equation" is interpreted as
# programmatically constructing the final string from its numeric components.
point1_ra_str = f"{ra1_h:02d} {ra1_m:02d} {ra1_s:02d}"
point1_dec_str = f"{dec1_decimal:.2f}"
point1_full_str = f"{point1_ra_str}, {point1_dec_str}"

point2_ra_str = f"{ra2_h:02d} {ra2_m:02d} {ra2_s:02d}"
point2_dec_str = f"{dec2_decimal:.2f}"
point2_full_str = f"{point2_ra_str}, {point2_dec_str}"

# Combine the two points into the final answer string, separated by a semicolon.
final_answer = f"{point1_full_str}; {point2_full_str}"

print(final_answer)