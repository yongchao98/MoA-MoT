import math

# Lengths of the sides
AB = math.sqrt(53)
BC = math.sqrt(65)
AC = 2

# Calculate the cosine of angle B
cos_B = (AB**2 + BC**2 - AC**2) / (2 * AB * BC)

# Calculate the angle in radians
angle_B_rad = math.acos(cos_B)

# Convert the angle to degrees
angle_B_deg = math.degrees(angle_B_rad)

# Round the result to 3 decimal places
angle_B_deg_rounded = round(angle_B_deg, 3)

print(angle_B_deg_rounded)