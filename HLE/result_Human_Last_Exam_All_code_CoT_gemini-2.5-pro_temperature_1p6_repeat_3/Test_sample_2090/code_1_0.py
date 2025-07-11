import math

# This script will determine and print the IAU boundary segment that the
# ecliptic most recently crossed when the Vernal Equinox moved into Pisces.
# The boundary is defined in the B1875.0 epoch.

# The boundary is a line segment at Right Ascension 0 hours.
# We define the coordinates of the two endpoints of this segment.

# Endpoint 1:
p1_ra_h = 0
p1_ra_m = 0
p1_ra_s = 0
p1_dec_deg = -6
p1_dec_min = 0

# Endpoint 2:
p2_ra_h = 0
p2_ra_m = 0
p2_ra_s = 0
p2_dec_deg = 5
p2_dec_min = 30

# Convert declination from degrees and arcminutes to decimal degrees.
# Formula: decimal_degrees = degrees + minutes / 60
p1_dec_decimal = p1_dec_deg + (p1_dec_min / 60.0)
p2_dec_decimal = p2_dec_deg + (p2_dec_min / 60.0)

# The user wants the numerically lowest point to appear first.
# We compare the decimal declinations to determine the order.
if p1_dec_decimal < p2_dec_decimal:
    lower_point_ra_h, lower_point_ra_m, lower_point_ra_s = p1_ra_h, p1_ra_m, p1_ra_s
    lower_point_dec = p1_dec_decimal
    
    higher_point_ra_h, higher_point_ra_m, higher_point_ra_s = p2_ra_h, p2_ra_m, p2_ra_s
    higher_point_dec = p2_dec_decimal
else:
    lower_point_ra_h, lower_point_ra_m, lower_point_ra_s = p2_ra_h, p2_ra_m, p2_ra_s
    lower_point_dec = p2_dec_decimal
    
    higher_point_ra_h, higher_point_ra_m, higher_point_ra_s = p1_ra_h, p1_ra_m, p1_ra_s
    higher_point_dec = p1_dec_decimal

# Format each point into the string "XX YY ZZ, AA.BB".
# We use f-strings for precise formatting, including leading zeros and two decimal places.
lower_point_str = (f"{lower_point_ra_h:02d} {lower_point_ra_m:02d} "
                   f"{lower_point_ra_s:02d}, {lower_point_dec:.2f}")

higher_point_str = (f"{higher_point_ra_h:02d} {higher_point_ra_m:02d} "
                    f"{higher_point_ra_s:02d}, {higher_point_dec:.2f}")

# Combine the two points with a semicolon as requested.
final_answer = f"{lower_point_str}; {higher_point_str}"

print(final_answer)
