# Parameters for the Goldberg Polyhedron determined from the image
# m is the number of steps along one axis between adjacent pentagons.
# n is the number of steps after a 60-degree turn.
# From the image, the path is straight, so n=0.
# There are 7 hexagons between two pentagons, so m = 7 + 1 = 8.
m = 8
n = 0

# The number of pentagonal faces (P) in a Goldberg polyhedron is always 12.
P = 12

# The triangulation number (T) is calculated as T = m^2 + mn + n^2.
T = m**2 + m * n + n**2

# The number of hexagonal faces (H) is calculated as H = 10 * (T - 1).
H = 10 * (T - 1)

# Print the final answer in the format m,n,H,P
print(f"{m},{n},{H},{P}")