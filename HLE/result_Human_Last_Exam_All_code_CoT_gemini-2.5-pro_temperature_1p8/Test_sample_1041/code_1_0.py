import math

# Define the parameters of the problem
R = 1000.0  # Radius of the outer disk
start_point = (0, 300)
target_point_1 = (0, 0)
target_point_2 = (2, 0)

# 1. Calculate the center of the target set
target_center_x = (target_point_1[0] + target_point_2[0]) / 2.0
target_center_y = (target_point_1[1] + target_point_2[1]) / 2.0
target_center = (target_center_x, target_center_y)

# 2. Calculate the distance 'd' from the start point to the target's center
d = math.sqrt((start_point[0] - target_center[0])**2 + (start_point[1] - target_center[1])**2)

# 3. Calculate the effective radius 'r_eff' (logarithmic capacity) of the target set
# For two points {z1, z2}, the capacity is |z1 - z2| / 2.
dist_between_targets = math.sqrt((target_point_1[0] - target_point_2[0])**2 + (target_point_1[1] - target_point_2[1])**2)
r_eff = dist_between_targets / 2.0

# 4. Apply the hitting probability formula: P = ln(R/d) / ln(R/r_eff)
numerator = math.log(R / d)
denominator = math.log(R / r_eff)
probability = numerator / denominator

# Output the explanation and the result
print("The probability is calculated using the formula from continuous potential theory:")
print("P = ln(R / d) / ln(R / r_eff)\n")
print("Where:")
print(f"R (outer radius) = {R}")
print(f"d (distance from start to target center) = sqrt(({start_point[0]}-{target_center[0]})^2 + ({start_point[1]}-{target_center[1]})^2) = {d:.4f}")
print(f"r_eff (effective radius of the target set) = |({target_point_2[0]},{target_point_2[1]}) - ({target_point_1[0]},{target_point_1[1]})| / 2 = {r_eff}")
print("\nCalculation:")
print(f"P = ln({R} / {d:.4f}) / ln({R} / {r_eff})")
print(f"P = {numerator:.4f} / {denominator:.4f}")
print(f"P = {probability:.4f}")
print("\nThe probability rounded to three significant digits is:")
print(f"{probability:.3f}")
