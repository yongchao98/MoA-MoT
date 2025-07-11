# Parameters m and n are determined by visual inspection of the polyhedron.
# By counting the steps along the hexagonal grid between two adjacent
# pentagons (centers of the colored spirals), we find a path
# of 4 steps in one direction, followed by a 60-degree turn
# and 2 steps in another direction.
# With the condition m > n, we have:
m = 4
n = 2

# The number of pentagonal faces in a Goldberg polyhedron is always 12.
P = 12

# The number of hexagonal faces is calculated using the formula:
# H = 10 * (T - 1), where T is the triangulation number.
# T = m^2 + m*n + n^2
T = m**2 + m*n + n**2
H = 10 * (T - 1)

# Print the final answer in the format m,n,H,P without spaces.
print(f"{m},{n},{H},{P}")