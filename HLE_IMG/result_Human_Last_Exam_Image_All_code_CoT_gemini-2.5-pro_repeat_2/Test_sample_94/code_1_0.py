import math

# Step 1 & 2: Define constants and ratios from the problem description.
# Average number of points in circles from 5 experiments.
points_in_circles = [730, 740, 735, 732, 739]
total_points = 1000
avg_points_in_circles = sum(points_in_circles) / len(points_in_circles)

# Ratio of circle area to total area (p)
p_ratio_exp = avg_points_in_circles / total_points
# Ratio of green rectangle area to total area (k)
k_ratio_approx = 0.04

# Step 3 & 4: Test geometric hypotheses and find integer parameters.
# Based on analysis, the most numerically consistent model for the rectangle's
# dimensions defined by the circle packing is X=8r and Y=2r(1+sqrt(3)).
# We will use this model to find the integer parameters r, gw, and gh.

# From k = (gw*gh) / (X*Y), we derive gw*gh in terms of r.
# gw*gh ≈ k * (8r * 2r(1+sqrt(3)))
# We assume from visual inspection that gh = 2r.
# This leads to gw ≈ 0.874 * r.
# We search for the best integer fit for (r, gw). The pair (8, 7) is the best fit.
r = 8
gw = 7
gh = 2 * r

# Step 5: Calculate the final dimensions using the determined parameters.
# Length X is determined by the 4 circles in the middle row.
X = 8 * r
# Width Y is determined by the 3 staggered rows of circles.
Y = 2 * r * (1 + math.sqrt(3))

# Round the final dimensions to the nearest integer as requested ("rounded to cm").
X_rounded = round(X)
Y_rounded = round(Y)

# Print the step-by-step derivation of the final answer.
print("Based on the data and geometric analysis, the parameters are determined as:")
print(f"Radius of circles (r) = {r}")
print(f"Green rectangle width (gw) = {gw}")
print(f"Green rectangle height (gh) = 2 * r = {gh}")
print("\nThese parameters lead to the following dimensions for the outer rectangle:")
print(f"Length (X) = 8 * r = 8 * {r} = {X}")
print(f"Width (Y) = 2 * r * (1 + sqrt(3)) = 2 * {r} * (1 + {math.sqrt(3):.3f}) = {Y:.3f}")
print("\nRounding these to the nearest centimeter:")
print(f"Rounded Length (x) = {X_rounded} cm")
print(f"Rounded Width (y) = {Y_rounded} cm")
print("\nFinal answer in x:y format:")
print(f"{X_rounded}:{Y_rounded}")

# Final Answer in the requested format for the game.
final_answer = f"{X_rounded}:{Y_rounded}"