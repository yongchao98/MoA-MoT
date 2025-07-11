# Step 1: Define the parameters m and n based on visual inspection of the polyhedron.
# By tracing the path between two adjacent pentagons (the centers of the colored spirals),
# we count 5 steps in a straight line, a 60-degree turn, and then 3 more steps.
# As per the hint, m > n.
m = 5
n = 3

# Step 2: The number of pentagonal faces (P) in a Goldberg polyhedron is always 12.
P = 12

# Step 3: Calculate the number of hexagonal faces (H).
# The formula for H is 10 * (T - 1), where T is the triangulation number.
# T = m^2 + m*n + n^2
T = m**2 + m * n + n**2
H = 10 * (T - 1)

# Step 4: Print the final answer in the specified format "m,n,H,P".
print(f"{m},{n},{H},{P}")