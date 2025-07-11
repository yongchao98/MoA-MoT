import math

# Step 1: Analyze Monte Carlo data
points_in_circles = [730, 740, 735, 732, 739]
total_points = 1000
avg_points_in_circles = sum(points_in_circles) / len(points_in_circles)
ratio_circles_to_total = avg_points_in_circles / total_points

# Step 2: Use the area relationship for the green rectangle
# A_g / A_total ≈ 0.04  => A_total / A_g ≈ 25
# We assume the exact relationship is 25 for integer-based geometry.
ratio_total_to_green = 25

# Step 3: Find the relationship between component dimensions (r, w_g, h_g)
# We have A_total = (10 * pi * r^2) / ratio_circles_to_total
# and A_total = 25 * A_g = 25 * w_g * h_g
# So, w_g * h_g / r^2 ≈ (10 * pi) / (25 * ratio_circles_to_total)
target_ratio_Ag_r_squared = (10 * math.pi) / (ratio_total_to_green * ratio_circles_to_total)

# Step 4: Search for the best integer solution for r and A_g
best_fit = {"r": 0, "A_g": 0, "error": float('inf')}
for r_candidate in range(1, 20):
    # From w_g*h_g ≈ target_ratio * r^2, find the nearest integer for A_g
    A_g_candidate = round(target_ratio_Ag_r_squared * (r_candidate**2))
    
    current_ratio = A_g_candidate / (r_candidate**2)
    error = abs(current_ratio - target_ratio_Ag_r_squared)
    
    if error < best_fit["error"]:
        best_fit["r"] = r_candidate
        best_fit["A_g"] = A_g_candidate
        best_fit["error"] = error

r = best_fit["r"]
A_g = best_fit["A_g"]

# Step 5: Based on the best-fit integers, calculate the definitive Total Area
# From the image, h_g seems to be 2*r. Let's find factors of A_g
# and see if this holds. For r=7, A_g=84. 84 = 6 * 14. h_g=14, 2*r=14. Match!
A_total = ratio_total_to_green * A_g

# Step 6: Determine L and W from geometry
# The width W appears to be 3 circle diameters stacked. W = 3 * (2*r) = 6*r
W = 6 * r
# The length L can be calculated from the total area and width
L = A_total / W

print("The final deduced values are:")
print(f"Circle radius (r) = {r}")
print(f"Green rectangle area (A_g) = {A_g}")
print(f"Total Area = {A_total}")
print(f"The length of the outer rectangle is {L:.0f} cm.")
print(f"The width of the outer rectangle is {W:.0f} cm.")
print("\nThe final answer is...")
print(f"{L:.0f}:{W:.0f}")
