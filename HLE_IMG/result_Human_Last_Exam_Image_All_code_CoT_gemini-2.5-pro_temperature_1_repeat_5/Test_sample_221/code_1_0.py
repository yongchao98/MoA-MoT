# Step 1: Define the parameters m and n based on visual inspection of the polyhedron.
# By tracing a path between two adjacent pentagons (centers of the colored spirals),
# we count 5 steps in one direction, a 60-degree turn, and 3 steps in the new direction.
# As per the hint, m > n.
m = 5
n = 3

# Step 2: Define the number of pentagonal faces, P.
# For any Goldberg polyhedron, the number of pentagons is always 12.
P = 12

# Step 3: Calculate the number of hexagonal faces, H, using the formula.
# The formula for H is 10 * (m^2 + m*n + n^2 - 1).
H = 10 * (m**2 + m*n + n**2 - 1)

# Step 4: Print the final answer in the format m,n,H,P without spaces.
print(f"{m},{n},{H},{P}")
