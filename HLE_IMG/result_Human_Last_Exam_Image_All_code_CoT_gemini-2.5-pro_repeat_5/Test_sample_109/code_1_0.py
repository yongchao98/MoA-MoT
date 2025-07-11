import math

# Step 1 & 2: Define initial parameters and geometric constraints
r_yellow = 1

# Step 3: Analyze the geometry of a yellow circle touching an edge and a large white circle.
# Let R be the radius of the large white circles.
# Let a yellow circle touch the top edge (y=H) and a large white circle in the top row.
# The y-center of the yellow circle is y_y = H - r = H - 1.
# The large white circles in the top row also touch the top edge, so their y-center is y_c = H - R.
# The vertical distance between their centers is d_y = y_y - y_c = (H - 1) - (H - R) = R - 1.
# The distance between their centers is R + r = R + 1.
# Let d_x be the horizontal distance between their centers.
# From the Pythagorean theorem: d_x^2 + d_y^2 = (R + r)^2
# d_x^2 + (R - 1)^2 = (R + 1)^2
# d_x^2 + R^2 - 2R + 1 = R^2 + 2R + 1
# d_x^2 = 4R
# For d_x to be an integer, R must be a perfect square. d_x = 2 * sqrt(R).

# Step 4: Analyze the row packing.
# Let g be the horizontal stagger and h_rows be the vertical distance between the centers of adjacent rows.
# For the circles in adjacent rows to touch, g^2 + h_rows^2 = (2R)^2.
# (g, h_rows, 2R) must be a Pythagorean triple.
# We test values of R that are perfect squares until we find one that can form a Pythagorean triple.
# R=1, 4, 9, 16 don't work.
# Let's test R = 25 (a perfect square, so d_x is an integer).
R = 25
d_x = 2 * math.sqrt(R) # d_x = 10
# We need to find integers g, h_rows such that g^2 + h_rows^2 = (2*25)^2 = 50^2 = 2500.
# The primitive Pythagorean triple (7, 24, 25) can be scaled by 2: (14, 48, 50).
# So, we can set g=14 and h_rows=48. This fits the visual representation.
g = 14
h_rows = 48

# Step 5: Determine the total height H.
# The middle row is at the vertical center of the image, H/2.
# The bottom row circles touch the bottom edge (y=0), so their center is at y=R.
# The vertical distance from the bottom row center to the middle row center is h_rows.
# So, H/2 - R = h_rows.
# H = 2 * (R + h_rows)
H = 2 * (R + h_rows)

# Step 6: Determine the width.
# The width of the circle pattern must be the same for all rows.
# Width of N circles of radius R placed side-by-side: N * 2R.
# Width of a C-Y-C arrangement (2 circles, 1 yellow separator): 2*R + 2*d_x.
# A C-C junction has width 2R. A C-Y-C arrangement replaces this with 2R + 2d_x, but the circles are separated by 2d_x instead of 2R.
# The width taken by N_T circles with k_T yellow separators is N_T*2R - k_T*(2R - 2*d_x).
# Let W_circles be the width of the circle pattern area.
# For the middle row (no yellow circles): W_circles = N_M * 2R
# For the top/bottom rows: W_circles = N_T*2R - k_T*(2*R - 2*d_x)
# Equating them: N_M * 2R = N_T*2R - k_T*(2*R - 2*d_x)
# N_M * 50 = N_T * 50 - k_T * (50 - 2*10) = N_T * 50 - k_T * 30
# 5 * N_M = 5 * N_T - 3 * k_T  => 5 * (N_T - N_M) = 3 * k_T
# For N_T and N_M to be integers, k_T must be a multiple of 5.
# Let's assume the simplest non-trivial case for the repeating pattern: k_T = 5.
# This implies 5 * (N_T - N_M) = 15 => N_T - N_M = 3.
# Let's use the simplest integers satisfying this and fitting the image: N_M=4, N_T=7.
# The width of the circle area is W_circles = N_M * 2*R = 4 * 2 * 25 = 200.
# The green bars are framing this circle area. The total target area consists of the circle pattern area and the green bar areas.
# To calculate the probability, we can consider a large representative section (a repeating unit) of the target.
# A repeating unit has a width of W_circles and height H.
W_circles = 4 * 2 * R
total_area_unit = W_circles * H

# Step 7: Calculate the total area of yellow circles in the repeating unit and the final expected value.
# Number of yellow circles in the top row of the unit is k_T = 5.
# Number of yellow circles in the bottom row is k_B = 5 (by symmetry).
# Number of yellow circles in the middle row is 0.
# The problem also shows yellow circles near the green bars, which are at the edges of the full target.
# The most robust calculation of probability comes from the density within the repeating circle pattern, as the green bar width is not specified and could be arbitrary.
# We consider the probability within the circle-filled area.
num_yellow_circles_in_unit = 5 + 5 # Top row + Bottom row
area_one_yellow_circle = math.pi * r_yellow**2
total_yellow_area_in_unit = num_yellow_circles_in_unit * area_one_yellow_circle

# The probability of hitting a yellow circle is the ratio of the yellow area to the total area.
probability = total_yellow_area_in_unit / total_area_unit

# The expected number of hits in 10000 shots.
shots = 10000
expected_hits = shots * probability

# Print the final equation with all numbers
print(f"Radius of large white circles (R): {R}")
print(f"Radius of small yellow circles (r): {r_yellow}")
print(f"Height of target (H): {H}")
print(f"Width of the repeating circle pattern (W_circles): {W_circles}")
print(f"Number of yellow circles in the repeating unit: {num_yellow_circles_in_unit}")
print("\nFinal Calculation:")
print(f"Expected Hits = {shots} * (Number of Yellow Circles * PI * r^2) / (W_circles * H)")
print(f"Expected Hits = {shots} * ({num_yellow_circles_in_unit} * PI * {r_yellow}^2) / ({W_circles} * {H})")
print(f"Expected Hits = {shots} * {num_yellow_circles_in_unit} * PI / {total_area_unit}")
print(f"Expected Hits = {10000 * num_yellow_circles_in_unit * math.pi / total_area_unit:.4f}")

# Final Answer
print(f"\nFinal Answer: {expected_hits}")