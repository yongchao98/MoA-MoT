import math

# Step 1: Use the area ratios to find the ratio of component dimensions.
# P_c = Area_circles / Area_total ≈ 0.7352
# P_g = Area_green / Area_total ≈ 0.04
# Area_circles = 10 * pi * r^2
# Area_green = w_g * h_g
# From the image, h_g = 2 * r
P_c = 0.7352
P_g = 0.04
ratio_of_areas = P_c / P_g

# (10 * pi * r^2) / (w_g * 2 * r) = ratio_of_areas
# (5 * pi * r) / w_g = ratio_of_areas
# w_g / r = (5 * pi) / ratio_of_areas
w_g_over_r = (5 * math.pi) / ratio_of_areas

# We know w_g and r are integers. We found that 6/7 is a good approximation.
r = 7
w_g = 6

# Step 2: Calculate all component dimensions.
h_g = 2 * r

# Step 3: Calculate total area from the green rectangle info.
area_green = w_g * h_g
area_total = area_green / P_g

# Step 4: Determine the outer rectangle dimensions.
# From the image, the total height y is 3 times the height of a component.
y = 3 * h_g
# From the total area, calculate x.
x = area_total / y

# Round to the nearest cm (they are already integers).
x_rounded = round(x)
y_rounded = round(y)

print(f"The analysis leads to the following dimensions:")
print(f"Circle radius (r): {r} cm")
print(f"Green rectangle width (w_g): {w_g} cm")
print(f"Green rectangle height (h_g): {h_g} cm")
print(f"Total Area: {area_total} cm^2")
print(f"Outer rectangle length (x): {x_rounded} cm")
print(f"Outer rectangle width (y): {y_rounded} cm")
print(f"\nThe final answer is in the format x:y")
print(f"{x_rounded}:{y_rounded}")
