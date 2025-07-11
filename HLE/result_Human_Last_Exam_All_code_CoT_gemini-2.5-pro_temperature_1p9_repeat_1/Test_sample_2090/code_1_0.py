import math

# The problem is to identify the specific IAU constellation boundary marker
# that the apparent position of the Sun at the Vernal Equinox crossed
# when it moved from Aries into Pisces due to precession.

# The IAU constellation boundaries are defined in the B1875.0 coordinate system.
# The boundary between Aries and Pisces is located at a Right Ascension of 0h 0m 0s.

# The Vernal Equinox is located on the celestial equator (Declination = 0).
# We need the two vertices of the boundary line that bracket the celestial equator.
# According to the official IAU definitions, these points are at
# Declination +5 degrees and -5 degrees.

# Point 1: The point with the numerically lower declination.
p1_ra_h = 0
p1_ra_m = 0
p1_ra_s = 0
p1_dec = -5.0

# Point 2: The point with the numerically higher declination.
p2_ra_h = 0
p2_ra_m = 0
p2_ra_s = 0
p2_dec = 5.0

# We will now format these two points into the required string format:
# XX YY ZZ, AA.BB; XX YY ZZ, AA.BB

# Format each point string
point1_str = f"{p1_ra_h:02d} {p1_ra_m:02d} {p1_ra_s:02d}, {p1_dec:.2f}"
point2_str = f"{p2_ra_h:02d} {p2_ra_m:02d} {p2_ra_s:02d}, {p2_dec:.2f}"

# Combine the two points with a semicolon as requested
# The problem asks to still output each number in the final equation.
# Here, the final "equation" is the line segment defined by the two points.
print(f"{point1_str}; {point2_str}")