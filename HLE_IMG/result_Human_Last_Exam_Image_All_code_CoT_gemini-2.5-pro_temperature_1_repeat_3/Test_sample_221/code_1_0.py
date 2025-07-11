# Step 1: Define the parameters m and n based on visual inspection of the polyhedron.
# The path between two pentagons involves 5 steps straight and 3 steps after a turn.
# Given m > n, we have m=5 and n=3.
m = 5
n = 3

# Step 2: The number of pentagonal faces (P) in a Goldberg icosahedron is always 12.
P = 12

# Step 3: Calculate the triangulation number (T).
T = m**2 + m * n + n**2

# Step 4: Calculate the number of hexagonal faces (H).
H = 10 * (T - 1)

# Step 5: Print the final answer in the format m,n,H,P without spaces.
print(f"{m},{n},{H},{P}")