import math

# Step 1: Define the model based on the visual information.
# L = 8r (from 4 tangent circles in the middle row spanning the length)
# W = k*r (where k is a factor to be determined)
# Area_total = L * W = 8k * r^2

# Step 2: Use the constraint about the green rectangles.
# We assume "around 4%" means exactly 4%, or 1/25.
# Area_green_rect / Area_total = 1/25
# (w_g * h_g) / (8k * r^2) = 1/25
# This gives a relationship: 25 * w_g * h_g = 8k * r^2
# Since w_g, h_g, and r are integers, this imposes a strong constraint on k.

# Step 3: Use the Monte Carlo data to estimate k.
# Average points in circles = (730 + 740 + 735 + 732 + 739) / 5 = 735.2
# Area_circles / Area_total ≈ 735.2 / 1000 = 0.7352
# (10 * pi * r^2) / (8k * r^2) ≈ 0.7352
# k ≈ (10 * pi) / (8 * 0.7352) ≈ 5.343

# Step 4: Find a rational k that fits the constraints.
# We need a rational number k ≈ 5.343 that satisfies 25 * w_g * h_g = 8k * r^2.
# Let's test the hypothesis k = 21/4 = 5.25, which is close to 5.343.
# Substituting k = 21/4 into the equation from Step 2:
# 25 * w_g * h_g = 8 * (21/4) * r^2
# 25 * w_g * h_g = 42 * r^2
# For w_g * h_g to be an integer, 42 * r^2 must be divisible by 25.
# Since 42 and 25 are coprime, r^2 must be a multiple of 25.
# This implies r must be a multiple of 5. The simplest non-trivial solution is r = 5.

# Step 5: Verify the solution and calculate the final dimensions.
# Let's confirm our chosen values.
r = 5
k = 21/4

# With r=5, the area of a green rectangle is:
# w_g_h_g = (42 * r**2) / 25 = (42 * 25) / 25 = 42.
# This is an integer, so it's a valid solution (e.g., dimensions 6cm x 7cm).

# The area ratio of circles with this model is:
# Ratio = (10 * pi) / (8 * k) = (10 * math.pi) / (8 * 21/4) = 5 * math.pi / 21 ≈ 0.748
# This is reasonably close to the experimental 0.7352, confirming our model is sound.

# Now, we calculate the dimensions of the outer rectangle.
L_val = 8 * r
W_val = k * r

print(f"The determined radius of the circles is r = {r} cm.")
print(f"The length of the rectangle is L = 8 * {r} = {L_val} cm.")
print(f"The width of the rectangle is W = (21/4) * {r} = {W_val} cm.")

# Step 6: Round and format the final answer as requested.
L_rounded = round(L_val)
W_rounded = round(W_val)

print(f"\nRounding the dimensions to the nearest cm:")
print(f"Length x = round({L_val}) = {L_rounded}")
print(f"Width y = round({W_val}) = {W_rounded}")

final_answer = f"{L_rounded}:{W_rounded}"
print(f"\nThe final answer in the format x:y is: {final_answer}")