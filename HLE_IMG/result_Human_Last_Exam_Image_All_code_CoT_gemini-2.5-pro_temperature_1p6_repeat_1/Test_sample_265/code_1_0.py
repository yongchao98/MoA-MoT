import math

# --- Step 1: Define Given Data and Geometry ---
P = 90 * 1000   # Force in N
tau_perm = 250    # Permissible shear stress in N/mm^2
b = 50          # Length of horizontal welds (W1, W2) in mm
d = 100         # Length of vertical weld (W3) in mm
e_dist = 150      # Eccentric distance from the end of the horizontal weld in mm

print("--- Problem Parameters ---")
print(f"Force (P): {P / 1000} kN")
print(f"Permissible Shear Stress (τ_perm): {tau_perm} N/mm^2")
print(f"Weld dimensions (b, d): {b} mm, {d} mm")
print("-" * 30 + "\n")

# --- Step 2: Calculate Weld Group Centroid ---
# Origin at the bottom-left corner of the C-shape
# Weld 1 (top): l1=50, centroid1=(25, 100)
# Weld 2 (bottom): l2=50, centroid2=(25, 0)
# Weld 3 (vertical): l3=100, centroid3=(0, 50)
L = 2 * b + d  # Total length of the weld
x_bar = (b * (b / 2) + b * (b / 2) + d * 0) / L
y_bar = (b * d + b * 0 + d * (d / 2)) / L

print("--- Weld Group Properties ---")
print(f"Total weld length (L): {L} mm")
print(f"Centroid location (x_bar, y_bar): ({x_bar:.2f} mm, {y_bar:.2f} mm)")

# --- Step 3: Calculate Moment (Torque) on Weld Group ---
# Eccentricity 'e' of the force P about the centroid
e = (b + e_dist) - x_bar
# Moment (Torque) T due to eccentricity
T = P * e
print(f"Eccentricity (e): {e:.2f} mm")
print(f"Moment/Torque (T): {T:.2f} N-mm")

# --- Step 4: Calculate Unit Polar Moment of Inertia (Ju) ---
# Moment of inertia for unit throat thickness (t=1mm)
# Weld 1 (top)
Ix1_u = (b * (d - y_bar)**2)
Iy1_u = (b**3 / 12) + (b * (b / 2 - x_bar)**2)
# Weld 2 (bottom)
Ix2_u = (b * (0 - y_bar)**2)
Iy2_u = (b**3 / 12) + (b * (b / 2 - x_bar)**2)
# Weld 3 (vertical)
Ix3_u = (d**3 / 12) + (d * (d / 2 - y_bar)**2)
Iy3_u = (d * (0 - x_bar)**2)

Ix_u = Ix1_u + Ix2_u + Ix3_u
Iy_u = Iy1_u + Iy2_u + Iy3_u
Ju = Ix_u + Iy_u
print(f"Unit Polar Moment of Inertia (Ju): {Ju:.2f} mm^3")
print("-" * 30 + "\n")

# --- Step 5: Calculate Stresses at the Critical Point ---
# The critical point is the corner farthest from the centroid. We will analyze
# the bottom-right corner (point D).
# Coordinates of point D relative to the centroid:
rx = b - x_bar
ry = 0 - y_bar

# Primary shear stress per unit thickness (acts downwards)
tau_py_unit = P / L

# Secondary shear stress components per unit thickness (from torque)
# Note: For a clockwise torque T, τ_sx = -T*ry/Ju and τ_sy = T*rx/Ju
tau_sx_unit = -T * ry / Ju
tau_sy_unit = T * rx / Ju

# Total stress components per unit thickness
tau_x_total_unit = tau_sx_unit
tau_y_total_unit = tau_sy_unit - tau_py_unit # Primary stress is downwards

# Resultant stress per unit thickness
tau_max_unit = math.sqrt(tau_x_total_unit**2 + tau_y_total_unit**2)

print("--- Stress Calculation (at critical point) ---")
print(f"Primary shear stress (τ_py) per unit thickness: {tau_py_unit:.2f} N/mm per mm")
print(f"Secondary shear stress (τ_sx) per unit thickness: {tau_sx_unit:.2f} N/mm per mm")
print(f"Secondary shear stress (τ_sy) per unit thickness: {tau_sy_unit:.2f} N/mm per mm")
print(f"Resultant max stress per unit thickness: {tau_max_unit:.2f} N/mm per mm")
print("-" * 30 + "\n")


# --- Step 6: Determine Required Weld Size ---
# The final equation is: Permissible Stress = Maximum Stress
# τ_perm = τ_max_unit / t
# Solve for required throat thickness 't'
t_req = tau_max_unit / tau_perm

# Convert throat thickness 't' to weld leg size 's' (t = s * 0.707)
s_req = t_req / 0.7071

print("--- Final Weld Size Calculation ---")
print("The required throat thickness (t) is calculated as:")
print(f"t = Maximum Stress per unit thickness / Permissible Stress")
print(f"t = {tau_max_unit:.2f} / {tau_perm:.1f} = {t_req:.2f} mm")
print("\nThe required weld leg size (s) is calculated from the throat thickness:")
print(f"s = t / 0.7071")
print(f"s = {t_req:.2f} / 0.7071 = {s_req:.2f} mm")

print(f"\n<<<The required size of the weld (leg) is {s_req:.2f} mm.>>>")
<<<14.10>>>