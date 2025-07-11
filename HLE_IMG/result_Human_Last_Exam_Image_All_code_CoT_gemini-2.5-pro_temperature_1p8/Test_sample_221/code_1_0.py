# Step 1: Define the parameters based on visual analysis of the polyhedron.
# The number of pentagons in a Goldberg polyhedron is always 12.
P = 12
# By inspecting the grid path between adjacent pentagons (spiral centers),
# we determine the parameters m and n. The path consists of 7 steps
# straight and 1 step after a turn.
m = 7
n = 1

# Step 2: Calculate the number of hexagonal faces (H).
# First, calculate the triangulation number, T.
T = m**2 + m*n + n**2

# Then, use T to find H.
H = 10 * (T - 1)

# Step 3: Print the results in the required format "m,n,H,P".
print(f"{m},{n},{H},{P}")