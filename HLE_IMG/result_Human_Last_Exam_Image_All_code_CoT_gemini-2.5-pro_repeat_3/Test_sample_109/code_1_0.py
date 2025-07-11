import math

# Step 1: Determine the radius of the white circles (R).
# Let R be the radius of a white circle and r_y be the radius of a yellow circle.
# We are given r_y = 1.
# In the top row, a yellow circle is between two white circles. All three are tangent to the top edge.
# The two white circles are also tangent to the yellow circle.
# Let's consider the horizontal distance between the centers of two adjacent circles (one white, one yellow)
# both tangent to a line. The formula is d_x = 2 * sqrt(R * r_y).
# For the two white circles to be tangent to each other through the yellow circle in the middle,
# the distance between their centers (2*R) must be equal to the sum of the horizontal distances
# from each to the yellow circle's center. So, 2*R = 2 * d_x = 2 * (2 * sqrt(R * r_y)).
# 2*R = 4 * sqrt(R * 1)
# R = 2 * sqrt(R)
# R^2 = 4 * R
# R * (R - 4) = 0
# Since R > 0, we have R = 4.
R = 4
r_y = 1

# Step 2: Determine the dimensions of the target (W, H) and the green bars (w_g, h_g).
# We assume all coordinates and measurements are integers.
# Let's analyze the tangencies on the right side of the target to find w_g and h_g.
# Let the origin (0,0) be the bottom-left corner of the target.
# The center of the rightmost white circle in the bottom row (W_b3) can be found by laying out the
# W-Y-W-Y-W sequence from the left edge. This places its center at (20, 4).
# This circle W_b3 is tangent to the bottom-right yellow circle (Y_r2).
# The center of Y_r2 is at (W-1, y_y2). The distance between their centers must be R + r_y = 5.
# ( (20) - (W-1) )^2 + (4 - y_y2)^2 = 5^2
# (21 - W)^2 + (4 - y_y2)^2 = 25
# The bottom-right green bar (G_b) has its bottom at y=0 and top at y=h_g. Y_r2 is tangent to its top.
# So, the center of Y_r2 is at y_y2 = h_g + r_y = h_g + 1.
# The right edge of W_b3 (at x=24) is tangent to the left edge of G_b (at x=W-w_g).
# So, 24 = W - w_g, which means W = 24 + w_g.
# Substitute W and y_y2 into the equation:
# (21 - (24 + w_g))^2 + (4 - (h_g + 1))^2 = 25
# (-3 - w_g)^2 + (3 - h_g)^2 = 25
# (w_g + 3)^2 + (h_g - 3)^2 = 25
# Since w_g and h_g are positive integers, we look for integer solutions. The pairs for (X^2, Y^2) that sum to 25 are (9, 16) and (16, 9).
# Case 1: (w_g+3)^2=9 -> w_g=0 or -6. (h_g-3)^2=16 -> h_g=7 or -1. Solution (w_g,h_g)=(0,7)
# Case 2: (w_g+3)^2=16 -> w_g=1 or -7. (h_g-3)^2=9 -> h_g=6 or 0. Solution (w_g,h_g)=(1,6)
# A green bar with width 0 is unlikely. The solution (w_g,h_g)=(1,6) is the most plausible.
w_g = 1
h_g = 6

# With w_g=1, we can find W.
W = 24 + w_g
# With h_g=6, we can find H.
# The middle row is at the vertical center of the image, H/2.
# The top-right yellow circle Y_r1 center's y-coordinate is H - h_g - r_y = H - 6 - 1 = H - 7.
# A key geometric deduction from the tangency of the rightmost middle-row white circle with Y_r1 shows their y-centers are identical.
# So, H/2 = H - 7, which gives H/2 = 7, so H = 14.
H = 14

# Step 3: Calculate the total area and the yellow area.
total_area = W * H

num_yellow_circles = 6
yellow_circle_radius = r_y
area_of_one_yellow_circle = math.pi * yellow_circle_radius**2
total_yellow_area = num_yellow_circles * area_of_one_yellow_circle

# Step 4: Calculate the expected number of hits for 10000 shots.
num_shots = 10000
probability_of_hit = total_yellow_area / total_area
expected_hits = num_shots * probability_of_hit

# Print the final calculation and result
print(f"The radius of a large white circle is R = {R} cm.")
print(f"The dimensions of a green bar are width = {w_g} cm and height = {h_g} cm.")
print(f"The dimensions of the target are width = {W} cm and height = {H} cm.")
print(f"Total area of the target = {W} * {H} = {total_area} cm^2.")
print(f"There are {num_yellow_circles} yellow circles, each with a radius of {yellow_circle_radius} cm.")
print(f"Total area of yellow circles = {num_yellow_circles} * pi * {yellow_circle_radius}^2 = {num_yellow_circles * math.pi:.2f} cm^2.")
print("\nThe expected number of hits in 10000 shots is calculated as:")
print(f"Expected Hits = {num_shots} * (Total Yellow Area / Total Target Area)")
print(f"Expected Hits = {num_shots} * ({num_yellow_circles} * pi * {yellow_circle_radius}^2) / ({W} * {H})")
print(f"Expected Hits = {expected_hits:.2f}")

# Final Answer
final_answer_value = 10000 * (6 * math.pi) / (25 * 14)
print(f"\nFinal Answer Value: {final_answer_value}")
# The problem asks for the answer in a specific format
# Let's provide the final numerical answer rounded to one decimal place as is common in such problems.
# <<<538.6>>>
final_answer_formatted = round(final_answer_value, 1)
