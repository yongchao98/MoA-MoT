import math

# Step 1: Explain the problem and the formula used.
# The number of grid squares 'k' that the triangle's perimeter passes through
# is the total number of times the perimeter crosses a vertical or horizontal grid line.
# This holds true provided the perimeter does not pass through any lattice points (grid intersections).
# The formula for k is the sum of the Manhattan distances between the integer-floored coordinates of the vertices.
# k = |floor(x2)-floor(x1)| + |floor(y2)-floor(y1)| +
#     |floor(x3)-floor(x2)| + |floor(y3)-floor(y2)| +
#     |floor(x1)-floor(x3)| + |floor(y1)-floor(y3)|
# where (x1, y1), (x2, y2), (x3, y3) are the coordinates of the triangle's vertices.

# Step 2: Determine the optimal placement of the triangle.
# To maximize k, we rotate the triangle. The optimal orientation is when the legs (sides of length 18)
# have slopes 1/2 and -2. This corresponds to using a direction vector based on the (1, 2) pair.
# We place the right-angle vertex V1 near the origin, at (epsilon, epsilon)
# to ensure no lattice points lie on the perimeter.

# Leg length
L = 18

# Define rotation using the (2,1) direction vector for one leg.
# cos(theta) = 2/sqrt(5), sin(theta) = 1/sqrt(5)
cos_theta = 2 / math.sqrt(5)
sin_theta = 1 / math.sqrt(5)

# Calculate the displacement vectors for the two legs from the right-angle vertex.
# u is the vector for the first leg
u_x = L * cos_theta
u_y = L * sin_theta
# w is the vector for the second leg, perpendicular to u
w_x = -L * sin_theta
w_y = L * cos_theta

# Step 3: Define vertex positions.
# We choose a small epsilon. A value of 0.01 is chosen as it's small enough to avoid issues
# but allows us to manipulate the floor function. The specific value helps maximize the result.
# w_x = -18/sqrt(5) which is approx -8.0498.
# floor(epsilon + w_x) = floor(0.01 - 8.0498) = floor(-8.0398) = -9. This choice is beneficial.
epsilon = 0.01

# Vertex 1 (right angle)
v1_x = epsilon
v1_y = epsilon

# Vertex 2 is V1 + u
v2_x = v1_x + u_x
v2_y = v1_y + u_y

# Vertex 3 is V1 + w
v3_x = v1_x + w_x
v3_y = v1_y + w_y

# Step 4: Calculate the integer-floored coordinates
f_v1_x = math.floor(v1_x)
f_v1_y = math.floor(v1_y)

f_v2_x = math.floor(v2_x)
f_v2_y = math.floor(v2_y)

f_v3_x = math.floor(v3_x)
f_v3_y = math.floor(v3_y)

# Step 5: Calculate k using the formula.
# Contribution from side V1V2
term1 = abs(f_v2_x - f_v1_x)
term2 = abs(f_v2_y - f_v1_y)
k_v1v2 = term1 + term2

# Contribution from side V2V3 (hypotenuse)
term3 = abs(f_v3_x - f_v2_x)
term4 = abs(f_v3_y - f_v2_y)
k_v2v3 = term3 + term4

# Contribution from side V3V1
term5 = abs(f_v1_x - f_v3_x)
term6 = abs(f_v1_y - f_v3_y)
k_v3v1 = term5 + term6

# Total k
k = k_v1v2 + k_v2v3 + k_v3v1

# Print the results step-by-step
print(f"Optimal placement uses a rotated triangle.")
print(f"Leg vectors (approx): u=({u_x:.2f}, {u_y:.2f}), w=({w_x:.2f}, {w_y:.2f})")
print(f"Triangle vertices (approximate):")
print(f"V1 = ({v1_x:.2f}, {v1_y:.2f})")
print(f"V2 = ({v2_x:.2f}, {v2_y:.2f})")
print(f"V3 = ({v3_x:.2f}, {v3_y:.2f})")
print("-" * 30)
print(f"Floored vertex coordinates:")
print(f"floor(V1) = ({f_v1_x}, {f_v1_y})")
print(f"floor(V2) = ({f_v2_x}, {f_v2_y})")
print(f"floor(V3) = ({f_v3_x}, {f_v3_y})")
print("-" * 30)
print("Calculating k as the sum of Manhattan distances of floored vertices:")
print(f"Contribution from leg 1 (V1-V2): |{f_v2_x} - {f_v1_x}| + |{f_v2_y} - {f_v1_y}| = {term1} + {term2} = {k_v1v2}")
print(f"Contribution from hypotenuse (V2-V3): |{f_v3_x} - {f_v2_x}| + |{f_v3_y} - {f_v2_y}| = {term3} + {term4} = {k_v2v3}")
print(f"Contribution from leg 2 (V3-V1): |{f_v1_x} - {f_v3_x}| + |{f_v1_y} - {f_v3_y}| = {term5} + {term6} = {k_v3v1}")
print("-" * 30)
print(f"Total number of squares k = {k_v1v2} + {k_v2v3} + {k_v3v1} = {k}")