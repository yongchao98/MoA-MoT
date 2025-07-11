import math

# Step 1: Define the coordinates of the boundary line segment from the IAU (Delporte/Roman) data.
# The boundary is a line of constant Right Ascension.
ra_h = 1
ra_m_decimal = 32.7

dec1_deg = 4
dec1_min = 3

dec2_deg = 10
dec2_min = 13

# Step 2: Convert Right Ascension to HH MM SS format.
# 32.7 minutes is 32 minutes and 0.7 * 60 = 42 seconds.
ra_m = int(ra_m_decimal)
ra_s = (ra_m_decimal - ra_m) * 60

# Format the RA string as XX YY ZZ
ra_str = f"{ra_h:02d} {ra_m:02d} {round(ra_s):02d}"

# Step 3: Convert Declinations to decimal degrees with two decimal places (AA.BB).
dec1_decimal = dec1_deg + dec1_min / 60.0
dec2_decimal = dec2_deg + dec2_min / 60.0

# Format the Declination strings
dec1_str = f"{dec1_decimal:.2f}"
dec2_str = f"{dec2_decimal:.2f}"

# Step 4: Assemble the final string for each point.
# The two points have the same RA but different declinations.
point1_str = f"{ra_str}, {dec1_str}"
point2_str = f"{ra_str}, {dec2_str}"

# The problem asks for the numerically lowest point first.
# Comparing the declinations, point 1 (4.05) is lower than point 2 (10.22).
final_answer = f"{point1_str}; {point2_str}"

print("The two points defining the boundary line segment are:")
print(final_answer)