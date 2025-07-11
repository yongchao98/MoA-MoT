import math

# Step 1: Define the triangle and the optimal orientation
L = 18.0
# The optimal orientation is found to be when tan(theta) = 1/2
theta = math.atan(0.5)
c = math.cos(theta)
s = math.sin(theta)

# Step 2: Define the coordinates of the vertices relative to the right-angle corner
# C is at (0,0), A and B are the other two vertices.
v_A = (L * c, L * s)
v_B = (-L * s, L * c)
v_C = (0.0, 0.0)

# We can place the triangle such that the right-angle corner C is at a point
# (epsilon, epsilon) for a very small epsilon > 0. This ensures no vertex or
# edge falls on a lattice point, and does not change the floor of non-integer coordinates.
x_coords = [v_C[0], v_A[0], v_B[0]]
y_coords = [v_C[1], v_A[1], v_B[1]]

print(f"Optimal orientation angle: theta = {math.degrees(theta):.2f} degrees")
print(f"Vertex A relative coordinates: ({v_A[0]:.3f}, {v_A[1]:.3f})")
print(f"Vertex B relative coordinates: ({v_B[0]:.3f}, {v_B[1]:.3f})")
print(f"Vertex C relative coordinates: ({v_C[0]:.3f}, {v_C[1]:.3f})")
print("-" * 20)

# Step 3: Calculate the floor of each coordinate
X_coords = [math.floor(x) for x in x_coords]
Y_coords = [math.floor(y) for y in y_coords]

print(f"Floored x-coordinates: {X_coords}")
print(f"Floored y-coordinates: {Y_coords}")
print("-" * 20)

# Step 4: Calculate the number of vertical (N_v) and horizontal (N_h) grid lines crossed
# The formula is 2 * (max_floor_coord - min_floor_coord)
N_v = 2 * (max(X_coords) - min(X_coords))
N_h = 2 * (max(Y_coords) - min(Y_coords))

print(f"Number of vertical crossings (N_v) = 2 * ({max(X_coords)} - ({min(X_coords)})) = {N_v}")
print(f"Number of horizontal crossings (N_h) = 2 * ({max(Y_coords)} - {min(Y_coords)}) = {N_h}")
print("-" * 20)

# Step 5: The total number of squares k is the sum of N_v and N_h
k = N_v + N_h
print(f"Total number of crossed squares k = N_v + N_h = {N_v} + {N_h} = {k}")
