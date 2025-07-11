# The coordinates for the boundary line vertices are based on the official
# IAU constellation boundaries defined by Eugène Delporte for the B1875.0 epoch.
# The boundary between Aries and Pisces is the line of Right Ascension = 0 hours.

# Right Ascension coordinates (same for both points)
ra_h = 0
ra_m = 0
ra_s = 0

# Declination for the two vertices of the boundary segment
# The segment spans from -10 degrees to +15 degrees, crossing the celestial equator.
dec_point1 = -10.00
dec_point2 = 15.00

# The problem asks for the numerically lowest point first.
# Comparing the declinations, -10.00 is lower than 15.00.
# So, we will list that point first.

# Format the output string as "XX YY ZZ, AA.BB; XX YY ZZ, AA.BB"
# Point 1 (lower declination)
point1_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {dec_point1:.2f}"

# Point 2 (higher declination)
point2_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}, {dec_point2:.2f}"

# Combine the two points with a semicolon
final_answer = f"{point1_str}; {point2_str}"

print(final_answer)