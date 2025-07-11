# The parameters m and n for a Goldberg polyhedron describe the path
# between the centers of two adjacent pentagonal faces on a hexagonal grid.
# By visual inspection of the image, we trace a path between a central
# pentagon (e.g., red) and an adjacent one (e.g., magenta).
# This path consists of 3 steps in one direction, a 60-degree turn,
# and 2 steps in another direction.
# Given the hint m > n, we assign:
m = 3
n = 2

# A Goldberg polyhedron always has 12 pentagonal faces.
P = 12

# The number of hexagonal faces (H) is calculated using the formula:
# H = 10 * (T - 1), where T = m^2 + mn + n^2.
# T is the triangulation number.
T = m**2 + m * n + n**2
H = 10 * (T - 1)

# Print the result in the specified format: m,n,H,P
print(f"{m},{n},{H},{P}")