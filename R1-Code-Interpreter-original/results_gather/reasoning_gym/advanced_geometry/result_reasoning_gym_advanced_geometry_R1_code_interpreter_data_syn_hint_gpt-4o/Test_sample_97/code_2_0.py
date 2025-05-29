import math

# Recalculate lengths of sides
AB = math.sqrt(241)
BC = math.sqrt(386)
AC = math.sqrt(17)

# Recalculate the angle at B using the law of cosines
cos_B = (AB**2 + BC**2 - AC**2) / (2 * AB * BC)
angle_B_rad = math.acos(cos_B)
angle_B_deg = math.degrees(angle_B_rad)

# Round the result to 3 decimal places
angle_B_deg_rounded = round(angle_B_deg, 3)

print(angle_B_deg_rounded)