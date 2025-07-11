# Step 1: Define the parameters m and n based on visual inspection of the polyhedron.
# By tracing the path between two adjacent pentagons (e.g., from the central red
# pentagon to the lower-right magenta one), we count m steps in one direction,
# a 60-degree turn, and n steps in the next.
# Visual counting reveals a path of 4 hexagons, a turn, and then 2 hexagons.
# Given the hint m > n, we set m=4 and n=2.
m = 4
n = 2

# Step 2: Define the number of pentagonal faces, P.
# For any Goldberg polyhedron, the number of pentagons is always 12.
P = 12

# Step 3: Calculate the number of hexagonal faces, H.
# First, calculate the triangulation number, T.
# The formula is T = m^2 + m*n + n^2.
T = m**2 + m*n + n**2

# Then, calculate the number of hexagons, H.
# The formula is H = 10 * (T - 1).
H = 10 * (T - 1)

# Step 4: Print the final answer in the format m,n,H,P.
print(f"{m},{n},{H},{P}")
