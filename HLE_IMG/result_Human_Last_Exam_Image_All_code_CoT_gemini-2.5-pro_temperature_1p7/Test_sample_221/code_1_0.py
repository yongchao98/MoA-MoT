# Parameters for the Goldberg Polyhedron determined from visual inspection.
# m represents the number of steps in the first straight line segment.
# n represents the number of steps after a 60-degree turn.
m = 4
n = 3

# The number of pentagonal faces in any Goldberg polyhedron is always 12.
P = 12

# The number of hexagonal faces is calculated using the formula:
# H = 10 * (m^2 + m*n + n^2 - 1)
H = 10 * (m**2 + m * n + n**2 - 1)

# The final output needs to be in the format m,n,H,P without spaces.
print(f"{m},{n},{H},{P}")