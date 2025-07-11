# Step 1: Define the parameters m and n based on visual inspection of the polyhedron.
# We count the steps between two adjacent pentagons.
# m is the number of steps in a straight line.
# n is the number of steps after a 60-degree turn.
m = 7
n = 4

# Step 2: The number of pentagonal faces (P) in a Goldberg polyhedron is always 12.
P = 12

# Step 3: Calculate the triangulation number (T).
# The formula for T is T = m^2 + m*n + n^2.
T = m**2 + m * n + n**2

# Step 4: Calculate the number of hexagonal faces (H).
# The formula for H is H = 10 * (T - 1).
H = 10 * (T - 1)

# Step 5: Print the final answer in the format m,n,H,P.
print(f"{m},{n},{H},{P}")