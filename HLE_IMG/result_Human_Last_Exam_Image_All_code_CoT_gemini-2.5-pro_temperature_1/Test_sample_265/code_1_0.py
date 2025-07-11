import math

# Step 1: Define given properties and weld geometry
P = 90000  # Force in N
tau_per = 250  # Permissible shear stress in N/mm^2

# Weld dimensions (in mm)
l1 = 50  # Length of top weld W1
l2 = 50  # Length of bottom weld W2
l3 = 100 # Length of vertical weld W3
e_dist = 150 # Eccentric distance from the end of the horizontal weld

# The welds are treated as lines for calculation.
# We define them with their length and centroid coordinates relative to an origin
# at the bottom-left corner of the vertical weld.
weld1 = {'length': l1, 'cx': l1 / 2, 'cy': l3}
weld2 = {'length': l2, 'cx': l2 / 2, 'cy': 0}
weld3 = {'length': l3, 'cx': 0, 'cy': l3 / 2}
welds = [weld1, weld2, weld3]

# --- Step 2: Calculate the Centroid of the Weld Group ---
A_total = l1 + l2 + l3
sum_Ax = weld1['length'] * weld1['cx'] + weld2['length'] * weld2['cx'] + weld3['length'] * weld3['cx']
sum_Ay = weld1['length'] * weld1['cy'] + weld2['length'] * weld2['cy'] + weld3['length'] * weld3['cy']

x_bar = sum_Ax / A_total
y_bar = sum_Ay / A_total

print("--- Weld Group Properties ---")
print(f"Total weld length (L): {A_total} mm")
print(f"Centroid location (x_bar, y_bar): ({x_bar:.2f} mm, {y_bar:.2f} mm)")

# Calculate eccentricity and Moment
force_x_location = l1 + e_dist
e = force_x_location - x_bar
M = P * e
print(f"Eccentricity (e): {force_x_location} mm - {x_bar:.2f} mm = {e:.2f} mm")
print(f"Moment (M = P * e): {P} N * {e:.2f} mm = {M:.2f} N-mm")
print("-" * 30)

# --- Step 3 & 5: Calculate Unit Moments of Inertia and Stresses ---
# Calculate Ixx and Iyy for the group about the centroid using Parallel Axis Theorem
# I = I_c + A*d^2, where I_c is moment of inertia about segment's own centroid.
# We calculate a "unit" moment of inertia assuming a throat thickness of 1 mm.
I_xx_u = 0
I_yy_u = 0

# Weld 1 (Top)
dy1 = weld1['cy'] - y_bar
dx1 = weld1['cx'] - x_bar
I_xx_u += (weld1['length'] * dy1**2) # I_c for horizontal line about x-axis is 0
I_yy_u += (weld1['length']**3 / 12) + (weld1['length'] * dx1**2)

# Weld 2 (Bottom)
dy2 = weld2['cy'] - y_bar
dx2 = weld2['cx'] - x_bar
I_xx_u += (weld2['length'] * dy2**2) # I_c for horizontal line about x-axis is 0
I_yy_u += (weld2['length']**3 / 12) + (weld2['length'] * dx2**2)

# Weld 3 (Vertical)
dy3 = weld3['cy'] - y_bar
dx3 = weld3['cx'] - x_bar
I_xx_u += (weld3['length']**3 / 12) + (weld3['length'] * dy3**2)
I_yy_u += (weld3['length'] * dx3**2) # I_c for vertical line about y-axis is 0

# Unit Polar Moment of Inertia
J_u = I_xx_u + I_yy_u

print("--- Moment of Inertia Calculation (for unit throat) ---")
print(f"Unit Moment of Inertia about X-axis (I_xx_u): {I_xx_u:.2f} mm^3")
print(f"Unit Moment of Inertia about Y-axis (I_yy_u): {I_yy_u:.2f} mm^3")
print(f"Unit Polar Moment of Inertia (J_u): {J_u:.2f} mm^3")
print("-" * 30)

# --- Step 4, 6: Find Maximum Resultant Stress ---
# The maximum stress will occur at the point on the weld furthest from the centroid.
# Let's check the top-right corner (point A) and bottom-right corner (point B).
# Point A: (50, 100), Point B: (50, 0)
point_A = {'x': l1, 'y': l3}

# Coordinates of Point A relative to the centroid
r_ax = point_A['x'] - x_bar
r_ay = point_A['y'] - y_bar

# Calculate stress components (as force per unit length, F')
# Primary shear force (acts vertically downwards)
F_prime_px = 0
F_prime_py = -P / A_total

# Secondary shear force due to torsion
F_prime_sx = - (M * r_ay) / J_u
F_prime_sy = (M * r_ax) / J_u

# Total shear force components
F_prime_total_x = F_prime_px + F_prime_sx
F_prime_total_y = F_prime_py + F_prime_sy

# Resultant shear force per unit length
F_prime_res = math.sqrt(F_prime_total_x**2 + F_prime_total_y**2)

print("--- Stress Calculation at Critical Point ---")
print(f"Primary Shear Force/Length (F'_p): {abs(F_prime_py):.2f} N/mm (downwards)")
print(f"Secondary Shear Force/Length Components (F'_sx, F'_sy): ({F_prime_sx:.2f}, {F_prime_sy:.2f}) N/mm")
print(f"Total Shear Force/Length Components (F'_total_x, F'_total_y): ({F_prime_total_x:.2f}, {F_prime_total_y:.2f}) N/mm")
print(f"Maximum Resultant Shear Force/Length (F'_res): {F_prime_res:.2f} N/mm")
print("-" * 30)

# --- Step 7: Determine the Weld Size ---
# The actual shear stress is tau = F_prime_res / t, where t is the throat thickness.
# We set tau = tau_per to find the required throat thickness.
# t = F'_res / tau_per
t = F_prime_res / tau_per

# For a fillet weld, throat thickness t = s * cos(45) = s * 0.7071
# So, the leg size s = t / 0.7071
s = t / 0.7071

print("--- Final Weld Size Calculation ---")
print("The required throat thickness (t) is found by the equation: t = F'_res / τ_per")
print(f"t = {F_prime_res:.2f} N/mm / {tau_per} N/mm^2 = {t:.4f} mm")
print("\nThe leg size (s) is found from the throat thickness: s = t / 0.7071")
print(f"s = {t:.4f} mm / 0.7071 = {s:.4f} mm")
print("\nTherefore, the required size of the weld (leg) is approximately {:.2f} mm.".format(s))

# Final answer in the required format
final_answer = s
#<<<14.11>>>