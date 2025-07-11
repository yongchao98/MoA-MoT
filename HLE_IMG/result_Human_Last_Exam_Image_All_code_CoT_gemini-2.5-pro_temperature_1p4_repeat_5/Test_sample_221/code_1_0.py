# Step 1: Define the parameters m and n based on visual inspection of the polyhedron.
# By tracing the path between two adjacent pentagons (e.g., the central one and the one below it),
# we find a path of 7 steps, a 60-degree turn, and 6 steps.
# As per the hint, m > n.
m = 7
n = 6

# Step 2: Define the number of pentagonal faces, P.
# For any Goldberg polyhedron, there are always 12 pentagonal faces.
P = 12

# Step 3: Calculate the number of hexagonal faces, H.
# First, calculate the triangulation number, T.
T = m**2 + m * n + n**2
# Then, use T to calculate H.
H = 10 * (T - 1)

# Step 4: Print the final answer in the format m,n,H,P.
print(f"{m},{n},{H},{P}")