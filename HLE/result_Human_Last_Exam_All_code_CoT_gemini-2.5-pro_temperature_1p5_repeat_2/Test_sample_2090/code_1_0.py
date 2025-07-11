import math

# The official IAU constellation boundaries are defined in B1875.0 coordinates.
# The boundary segment where the ecliptic crossed from Aries into Pisces
# is a line of constant Right Ascension.

# Right Ascension for the boundary line
ra_h = 1
ra_m = 59
ra_s = 48

# Declinations for the two endpoints of the line segment
# Endpoint 1: +11 degrees 46 arcminutes
# Endpoint 2: +15 degrees 16 arcminutes

# Convert declination to decimal degrees, rounded to two places
dec1_deg = 11
dec1_min = 46
dec1_decimal = round(dec1_deg + dec1_min / 60.0, 2)

dec2_deg = 15
dec2_min = 16
dec2_decimal = round(dec2_deg + dec2_min / 60.0, 2)

# The point with the lower declination value is listed first.
point1_ra_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}"
point1_dec_str = f"{dec1_decimal:.2f}"
point1_str = f"{point1_ra_str}, {point1_dec_str}"

point2_ra_str = f"{ra_h:02d} {ra_m:02d} {ra_s:02d}"
point2_dec_str = f"{dec2_decimal:.2f}"
point2_str = f"{point2_ra_str}, {point2_dec_str}"

# Print the final formatted string
# Each number is outputted to form the final string
print(f"{point1_ra_str}, {point1_dec_str}; {point2_ra_str}, {point2_dec_str}")
