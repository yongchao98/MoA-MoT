import math

# --- Given Parameters ---
P = 90000.0  # Force in N (90 kN)
tau_perm = 250.0  # Permissible shear stress in N/mm^2
l1 = 50.0  # Length of top horizontal weld (W1) in mm
l2 = 50.0  # Length of bottom horizontal weld (W2) in mm
l3 = 100.0 # Length of vertical weld (W3) in mm
e_dist = 150.0 # Distance from weld end to force application in mm

print("--- Step 1: Weld Group Properties and Centroid Calculation ---")
# Define weld segments with their length and centroid coordinates (origin at bottom-left corner of W3)
weld1 = {'length': l1, 'cx': 25.0, 'cy': 100.0}
weld2 = {'length': l2, 'cx': 25.0, 'cy': 0.0}
weld3 = {'length': l3, 'cx': 0.0, 'cy': 50.0}
welds = [weld1, weld2, weld3]

# Calculate total length and centroid of the weld group
L_total = sum(w['length'] for w in welds)
x_bar = sum(w['length'] * w['cx'] for w in welds) / L_total
y_bar = sum(w['length'] * w['cy'] for w in welds) / L_total

print(f"Total weld length (L_total): {L_total:.2f} mm")
print(f"Centroid of weld group (G): ({x_bar:.2f} mm, {y_bar:.2f} mm)")

print("\n--- Step 2: Load Calculation (Force and Moment) ---")
# Eccentricity (distance from centroid to force)
force_x_location = l1 + e_dist
e = force_x_location - x_bar
# Torque (Moment) about the centroid
T = P * e
print(f"Eccentricity (e): {force_x_location:.2f} - {x_bar:.2f} = {e:.2f} mm")
print(f"Moment (Torque, T): {P:.0f} N * {e:.2f} mm = {T:.2f} N-mm")

print("\n--- Step 3: Unit Polar Moment of Inertia (J_u) Calculation ---")
# Calculate unit moments of inertia (I_u) for each segment about the group centroid G
# and sum them up. Using parallel axis theorem: I_G = I_c + A*d^2
# For a line, Area A is replaced by length L.
# W1 (horizontal):
Ix_u1 = (l1 * (weld1['cy'] - y_bar)**2)
Iy_u1 = (l1**3 / 12) + (l1 * (weld1['cx'] - x_bar)**2)
# W2 (horizontal):
Ix_u2 = (l2 * (weld2['cy'] - y_bar)**2)
Iy_u2 = (l2**3 / 12) + (l2 * (weld2['cx'] - x_bar)**2)
# W3 (vertical):
Ix_u3 = (l3**3 / 12) + (l3 * (weld3['cy'] - y_bar)**2)
Iy_u3 = (l3 * (weld3['cx'] - x_bar)**2)

# Total unit moments of inertia
Ix_u_total = Ix_u1 + Ix_u2 + Ix_u3
Iy_u_total = Iy_u1 + Iy_u2 + Iy_u3
# Total unit polar moment of inertia
J_u = Ix_u_total + Iy_u_total
print(f"Unit Moment of Inertia about x-axis (Ix_u): {Ix_u_total:.2f} mm^3")
print(f"Unit Moment of Inertia about y-axis (Iy_u): {Iy_u_total:.2f} mm^3")
print(f"Unit Polar Moment of Inertia (J_u): {J_u:.2f} mm^3")

print("\n--- Step 4: Stress Calculation at Critical Point ---")
# The critical point is the top-right corner A(50, 100)
crit_point_x = 50.0
crit_point_y = 100.0
# Coordinates of critical point relative to centroid G
rx = crit_point_x - x_bar
ry = crit_point_y - y_bar
print(f"Critical point A located at ({crit_point_x:.1f}, {crit_point_y:.1f}).")
print(f"Position vector from G to A: r = ({rx:.2f}, {ry:.2f}) mm")

# Calculate forces per unit length (stresses are these values divided by throat thickness t)
# Primary shear force (acts vertically downwards)
f_primary_y = P / L_total
# Secondary (torsional) shear forces
# Note: Torque T is counter-clockwise.
# f_secondary_x = -T*ry/J_u (points left for positive ry)
# f_secondary_y = +T*rx/J_u (points up for positive rx)
f_secondary_x = T * ry / J_u
f_secondary_y = T * rx / J_u

print(f"Primary shear force per unit length (f'_y): {P:.0f} N / {L_total:.2f} mm = {f_primary_y:.2f} N/mm (downward)")
print(f"Secondary shear (x-comp, f''_x): ({T:.0f} * {ry:.2f}) / {J_u:.2f} = {f_secondary_x:.2f} N/mm (leftward)")
print(f"Secondary shear (y-comp, f''_y): ({T:.0f} * {rx:.2f}) / {J_u:.2f} = {f_secondary_y:.2f} N/mm (upward)")

# Vector sum of forces per unit length
f_total_x = -f_secondary_x # acting left
f_total_y = -f_primary_y + f_secondary_y # primary is down, secondary is up

f_resultant = math.sqrt(f_total_x**2 + f_total_y**2)

print("\n--- Step 5: Final Weld Size Calculation ---")
print("Equation for resultant force per unit length (f_res):")
print(f"f_res = sqrt( (f''_y - f'_y)^2 + (-f''_x)^2 )")
print(f"f_res = sqrt( ({f_secondary_y:.2f} - {f_primary_y:.2f})^2 + ({-f_secondary_x:.2f})^2 )")
print(f"f_res = sqrt( ({f_total_y:.2f})^2 + ({f_total_x:.2f})^2 ) = {f_resultant:.2f} N/mm")

# Required throat thickness (t)
# tau_perm = f_resultant / t
t_throat = f_resultant / tau_perm
print("\nEquation for required throat thickness (t):")
print(f"t = f_res / τ_perm = {f_resultant:.2f} N/mm / {tau_perm:.2f} N/mm^2 = {t_throat:.3f} mm")

# Required leg size (s) for a fillet weld
s_leg = t_throat / math.cos(math.radians(45)) # or t / 0.707
print("\nEquation for required weld leg size (s):")
print(f"s = t / cos(45°) = {t_throat:.3f} mm / {math.cos(math.radians(45)):.3f} = {s_leg:.2f} mm")

print(f"\nThe required weld leg size is {s_leg:.2f} mm.")
