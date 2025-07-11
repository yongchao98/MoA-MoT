import math

# Plan:
# 1. Define the coordinates for the two endpoints of the boundary line segment.
#    This segment is at a constant Right Ascension of 0 hours.
#    The data is from the official B1875.0 epoch IAU constellation boundary definitions.
# 2. Convert the declination values from degrees/arcminutes to decimal degrees.
# 3. Format the two points according to the specified format "XX YY ZZ, AA.BB".
# 4. Print the final string, ensuring the point with the lower declination value is first.

# Point 1 coordinates
ra_h_1, ra_m_1, ra_s_1 = 0, 0, 0
dec_deg_1, dec_min_1 = -5, 45

# Point 2 coordinates
ra_h_2, ra_m_2, ra_s_2 = 0, 0, 0
dec_deg_2, dec_min_2 = 5, 45

# Convert declinations to decimal degrees
# For negative declination, both parts are negative: -5 deg 45 min = -(5 + 45/60)
dec_decimal_1 = float(dec_deg_1) - float(dec_min_1) / 60.0
# For positive declination: 5 deg 45 min = 5 + 45/60
dec_decimal_2 = float(dec_deg_2) + float(dec_min_2) / 60.0

# Ensure the numerically lowest point is first by checking declinations
# In this case, dec_decimal_1 (-5.75) is lower than dec_decimal_2 (5.75)
point_low = {
    "ra_h": ra_h_1, "ra_m": ra_m_1, "ra_s": ra_s_1, "dec": dec_decimal_1
}
point_high = {
    "ra_h": ra_h_2, "ra_m": ra_m_2, "ra_s": ra_s_2, "dec": dec_decimal_2
}

# Format the output strings for each number and symbol
# This satisfies the requirement "output each number in the final equation!"
p1_ra_h_str = f"{point_low['ra_h']:02d}"
p1_ra_m_str = f"{point_low['ra_m']:02d}"
p1_ra_s_str = f"{point_low['ra_s']:02d}"
p1_dec_str = f"{point_low['dec']:.2f}"

p2_ra_h_str = f"{point_high['ra_h']:02d}"
p2_ra_m_str = f"{point_high['ra_m']:02d}"
p2_ra_s_str = f"{point_high['ra_s']:02d}"
p2_dec_str = f"{point_high['dec']:.2f}"

separator_coords = ", "
separator_points = "; "

# Print the final result by assembling the formatted parts
print(
    p1_ra_h_str + " " +
    p1_ra_m_str + " " +
    p1_ra_s_str +
    separator_coords +
    p1_dec_str +
    separator_points +
    p2_ra_h_str + " " +
    p2_ra_m_str + " " +
    p2_ra_s_str +
    separator_coords +
    p2_dec_str
)