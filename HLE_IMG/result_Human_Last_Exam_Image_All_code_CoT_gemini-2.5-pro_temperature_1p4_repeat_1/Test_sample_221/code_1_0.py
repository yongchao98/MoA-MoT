# Parameters for the Goldberg Polyhedron determined from visual analysis
m = 4
n = 1

# The number of pentagonal faces (P) in any Goldberg polyhedron is always 12.
P = 12

# The number of hexagonal faces (H) is calculated using the formula:
# H = 10 * (T - 1), where T is the triangulation number.
# T = m^2 + m*n + n^2

# Calculate the triangulation number T
T = m**2 + m * n + n**2

# Calculate the number of hexagonal faces H
H = 10 * (T - 1)

# Print the final answer in the format m,n,H,P without spaces
print(f"{m},{n},{H},{P}")