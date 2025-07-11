# A python program to calculate the parameters of the Goldberg polyhedron.

# Step 1: Define m and n based on visual inspection of the image.
# The path between adjacent pentagons consists of 5 steps straight and 2 steps after a turn.
# Given the hint m > n, we have m=5 and n=2.
m = 5
n = 2

# Step 2: Define the number of pentagonal faces, P.
# For any Goldberg polyhedron, P is always 12.
P = 12

# Step 3: Calculate the number of hexagonal faces, H.
# The formula is H = 10 * (T - 1), where T = m^2 + m*n + n^2.
# First, calculate the triangulation number, T.
T = m**2 + m * n + n**2
# Then, calculate H.
H = 10 * (T - 1)

# Step 4: Print the final answer in the format m,n,H,P
print(f"{m},{n},{H},{P}")