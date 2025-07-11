import math

# --- Given Data ---
P = 90 * 1000  # Force in N
tau_perm = 250  # Permissible shear stress in N/mm^2
l_w1 = 50  # Length of weld W1 in mm
l_w2 = 50  # Length of weld W2 in mm
l_w3 = 100 # Length of weld W3 in mm
d_lever = 150 # Lever arm from weld end to force in mm

print("--- Step 1: Analyze Weld Geometry ---")
# Coordinate system origin at the bottom-left corner of the weld group.
# Weld W1: horizontal, from (0, 100) to (50, 100)
# Weld W2: horizontal, from (0, 0) to (50, 0)
# Weld W3: vertical, from (0, 0) to (0, 100)
L_total = l_w1 + l_w2 + l_w3
print(f"Total weld length (L): {L_total} mm")
print("-" * 30)

print("--- Step 2: Locate Centroid (G) ---")
# Centroids of individual welds
x1_c, y1_c = l_w1 / 2, 100
x2_c, y2_c = l_w2 / 2, 0
x3_c, y3_c = 0, l_w3 / 2

# Centroid of the entire weld group
x_bar = (l_w1 * x1_c + l_w2 * x2_c + l_w3 * x3_c) / L_total
y_bar = (l_w1 * y1_c + l_w2 * y2_c + l_w3 * y3_c) / L_total
print(f"Centroid coordinates (x_bar, y_bar): ({x_bar:.2f} mm, {y_bar:.2f} mm)")
print("-" * 30)

print("--- Step 3: Calculate Torque ---")
# Eccentricity is the distance from the centroid to the line of action of the force
force_x_pos = l_w1 + d_lever
e = force_x_pos - x_bar
T = P * e
print(f"Force (P): {P} N")
print(f"Eccentricity (e): {force_x_pos} - {x_bar:.2f} = {e:.2f} mm")
print(f"Torque (T): {P} N * {e:.2f} mm = {T:.2f} N-mm")
print("-" * 30)

print("--- Step 4: Calculate Polar Moment of Inertia (per unit throat) ---")
# Using Parallel Axis Theorem: I = I_c + A*d^2
# I_xx (about the group's centroidal x-axis)
I_c_xx1 = 0  # For horizontal line
d_y1 = y1_c - y_bar
I_xx1 = I_c_xx1 + l_w1 * d_y1**2

I_c_xx2 = 0  # For horizontal line
d_y2 = y2_c - y_bar
I_xx2 = I_c_xx2 + l_w2 * d_y2**2

I_c_xx3 = l_w3**3 / 12 # For vertical line
d_y3 = y3_c - y_bar
I_xx3 = I_c_xx3 + l_w3 * d_y3**2

I_xx_u = I_xx1 + I_xx2 + I_xx3
print(f"Moment of inertia I_xx_u = {I_xx1:.2f} + {I_xx2:.2f} + {I_xx3:.2f} = {I_xx_u:.2f} mm^3")

# I_yy (about the group's centroidal y-axis)
I_c_yy1 = l_w1**3 / 12 # For horizontal line
d_x1 = x1_c - x_bar
I_yy1 = I_c_yy1 + l_w1 * d_x1**2

I_c_yy2 = l_w2**3 / 12 # For horizontal line
d_x2 = x2_c - x_bar
I_yy2 = I_c_yy2 + l_w2 * d_x2**2

I_c_yy3 = 0 # For vertical line
d_x3 = x3_c - x_bar
I_yy3 = I_c_yy3 + l_w3 * d_x3**2

I_yy_u = I_yy1 + I_yy2 + I_yy3
print(f"Moment of inertia I_yy_u = {I_yy1:.2f} + {I_yy2:.2f} + {I_yy3:.2f} = {I_yy_u:.2f} mm^3")

J_u = I_xx_u + I_yy_u
print(f"Polar moment of inertia J_u = {I_xx_u:.2f} + {I_yy_u:.2f} = {J_u:.2f} mm^3")
print("-" * 30)

print("--- Step 5 & 6: Find Maximum Resultant Stress ---")
# Critical points are the corners/ends of the weld.
# Point A: (50, 100), Point D: (50, 0) are farthest and most likely to be critical.
critical_points = {"A": (50, 100), "B": (0, 100), "C": (0, 0), "D": (50, 0)}
max_f_res = 0
critical_point_name = ""

# Primary shear force per unit length (acts in -y direction)
f_p_y = -P / L_total

for name, (x, y) in critical_points.items():
    # Secondary shear force components per unit length
    f_s_x = T * (y - y_bar) / J_u
    f_s_y = -T * (x - x_bar) / J_u

    # Total resultant force components per unit length
    f_res_x = f_s_x
    f_res_y = f_p_y + f_s_y
    
    # Resultant force per unit length
    f_res = math.sqrt(f_res_x**2 + f_res_y**2)
    
    print(f"At Point {name}({x},{y}): Resultant force/length = {f_res:.2f} N/mm")

    if f_res > max_f_res:
        max_f_res = f_res
        critical_point_name = name

print(f"\nMaximum resultant force per unit length (f_max) occurs at Point {critical_point_name} and is {max_f_res:.2f} N/mm")
print("-" * 30)

print("--- Step 7: Calculate Required Weld Size ---")
# Required throat thickness (t)
# τ_max = f_max / t  =>  t = f_max / τ_perm
t_req = max_f_res / tau_perm
print(f"Required throat thickness (t) = {max_f_res:.2f} N/mm / {tau_perm} N/mm^2 = {t_req:.2f} mm")

# Required leg size (s)
# t = s * sin(45) => s = t / sin(45)
s_req = t_req / math.sin(math.radians(45))
print(f"Required leg size (s) = {t_req:.2f} mm / sin(45) = {s_req:.2f} mm")
print("-" * 30)
print(f"Final Answer: The required weld size (leg) is {s_req:.2f} mm.")
