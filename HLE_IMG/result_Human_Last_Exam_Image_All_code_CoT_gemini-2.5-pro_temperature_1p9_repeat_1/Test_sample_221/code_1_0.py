import math

# Step 1: Determine the parameters from visual inspection of the polyhedron.
# By tracing the path between the centers of two adjacent colored spirals (pentagons),
# we move 4 steps in one direction, turn, and move 2 steps in the new direction.
# Given the hint m > n, we have:
m = 4
n = 2

# Step 2: The number of pentagonal faces in a Goldberg polyhedron is always 12.
P = 12

# Step 3: Calculate the number of hexagonal faces (H) using the standard formula.
# The triangulation number T is m^2 + m*n + n^2.
T = m**2 + m * n + n**2
# The number of hexagonal faces H is 10 * T - 10.
H = 10 * T - 10

# Step 4: Print the final answer in the required format "m,n,H,P".
# The final equation components are:
# m = 4
# n = 2
# H = 10 * (4^2 + 4*2 + 2^2) - 10 = 10 * (16 + 8 + 4) - 10 = 10 * 28 - 10 = 270
# P = 12
print(f"{m},{n},{H},{P}")
