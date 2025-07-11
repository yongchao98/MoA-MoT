import math

# Step 1: Define parameters m and n based on visual analysis of the polyhedron.
# The path between adjacent pentagons involves m steps, a 60-degree turn, and n steps.
# Visual inspection of the image reveals a path of 4 steps and 2 steps.
# Given the hint m > n, we have:
m = 4
n = 2

# Step 2: Define the number of pentagonal faces, P.
# This is a constant for all Goldberg polyhedra.
P = 12

# Step 3: Calculate the number of hexagonal faces, H.
# The formula is H = 10 * (T - 1), where T = m^2 + m*n + n^2.
T = m**2 + m*n + n**2
H = 10 * (T - 1)

# Step 4: Print the results in the required format "m,n,H,P".
print(f"{m},{n},{H},{P}")