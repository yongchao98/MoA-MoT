import math

# 1. Input parameters from the problem description
P = 90 * 1000  # Force in N (90 kN)
tau_perm = 250   # Permissible shear stress in N/mm^2
l_hor = 50       # Length of horizontal welds (W1, W2) in mm
l_ver = 100      # Length of vertical weld (W3) in mm
d = 150          # Additional distance for eccentricity in mm

print("--- Given Parameters ---")
print(f"Force, P = {P / 1000} kN = {P} N")
print(f"Permissible shear stress, τ_perm = {tau_perm} N/mm²")
print(f"Weld dimensions: l1 = {l_hor} mm, l2 = {l_hor} mm, l3 = {l_ver} mm")
print(f"Eccentricity distance, d = {d} mm\n")

# --- Step 1: Find the Centroid of the Weld Group ---
# We assume unit throat thickness (t=1mm) for these calculations.
# The "area" A is therefore the length of the weld.
# We set the origin at the bottom-left corner of the weld group.
# Weld W1: from (0, 100) to (50, 100). Centroid (25, 100), Area A1=50
# Weld W2: from (0, 0) to (50, 0). Centroid (25, 0), Area A2=50
# Weld W3: from (0, 0) to (0, 100). Centroid (0, 50), Area A3=100
A1, x1, y1 = l_hor, l_hor / 2, l_ver
A2, x2, y2 = l_hor, l_hor / 2, 0
A3, x3, y3 = l_ver, 0, l_ver / 2

A_total = A1 + A2 + A3
x_bar = (A1 * x1 + A2 * x2 + A3 * x3) / A_total
y_bar = (A1 * y1 + A2 * y2 + A3 * y3) / A_total

print("--- Step 1: Weld Group Centroid (G) ---")
print(f"Total weld length (acting as unit area), A_total = {A_total:.2f} mm")
print(f"Centroid coordinates, G = (x̄, ȳ) = ({x_bar:.2f} mm, {y_bar:.2f} mm)\n")

# --- Step 2: Calculate Polar Moment of Inertia (per unit throat) ---
# Using the parallel axis theorem: I_G = I_c + A*d^2, with respect to the weld group centroid G(x_bar, y_bar)
# For weld W1 (horizontal)
Ix1 = (A1 * (y1 - y_bar)**2)
Iy1 = (l_hor**3 / 12) + (A1 * (x1 - x_bar)**2)
# For weld W2 (horizontal)
Ix2 = (A2 * (y2 - y_bar)**2)
Iy2 = (l_hor**3 / 12) + (A2 * (x2 - x_bar)**2)
# For weld W3 (vertical)
Ix3 = (l_ver**3 / 12) + (A3 * (y3 - y_bar)**2)
Iy3 = (A3 * (x3 - x_bar)**2)

I_x_total = Ix1 + Ix2 + Ix3
I_y_total = Iy1 + Iy2 + Iy3

# Polar moment of inertia Ju for unit throat thickness
J_u = I_x_total + I_y_total

print("--- Step 2: Polar Moment of Inertia (per unit throat thickness) ---")
print(f"Moment of inertia about x-axis, Ix = {I_x_total:.2f} mm³")
print(f"Moment of inertia about y-axis, Iy = {I_y_total:.2f} mm³")
print(f"Polar moment of inertia, Ju = Ix + Iy = {J_u:.2f} mm³\n")

# --- Step 3: Calculate Stresses ---
# Eccentricity is the horizontal distance from the force P to the centroid G
eccentricity = (l_hor + d) - x_bar
M = P * eccentricity

# Primary shear stress (per unit throat), τ' - acts vertically downwards
tau_p_u = P / A_total

# Secondary shear stress (per unit throat), τ'' - due to moment M
# Maximum stress occurs at the point farthest from the centroid. Let's analyze the top-right corner (50, 100).
# Coordinates of critical point relative to the centroid G
x_r = l_hor - x_bar
y_r = l_ver - y_bar

# Secondary stress components at the critical point
tau_s_x_u = -(M * y_r) / J_u
tau_s_y_u = (M * x_r) / J_u

print("--- Step 3: Stress Calculations (per unit throat thickness) ---")
print(f"Eccentricity, e = {eccentricity:.2f} mm")
print(f"Bending Moment, M = P * e = {M:.2f} N-mm")
print(f"Primary Shear Stress (downward), τ' = P / A_total = {tau_p_u:.2f} N/mm")
print("\nAt critical point (top right corner):")
print(f"Secondary Shear Stress x-component, τ''x = -M*y_r/Ju = {tau_s_x_u:.2f} N/mm")
print(f"Secondary Shear Stress y-component, τ''y = M*x_r/Ju = {tau_s_y_u:.2f} N/mm")

# --- Step 4: Calculate Maximum Resultant Stress ---
# Combine stress components vectorially
tau_total_x = tau_s_x_u
tau_total_y = tau_s_y_u - tau_p_u  # τ' is downward, hence subtracted from τ''y

# Resultant stress (per unit throat)
tau_res_u = math.sqrt(tau_total_x**2 + tau_total_y**2)

print("\n--- Step 4: Maximum Resultant Stress ---")
print(f"Total Shear x-component, τ_x = {tau_total_x:.2f} N/mm")
print(f"Total Shear y-component, τ_y = {tau_s_y_u:.2f} - {tau_p_u:.2f} = {tau_total_y:.2f} N/mm")
print(f"Resultant Shear (per unit throat), τ_res_u = sqrt(τ_x² + τ_y²) = {tau_res_u:.2f} N/mm\n")

# --- Step 5: Determine Required Weld Size ---
# Required throat thickness t: τ_actual = τ_res_u / t <= τ_perm
t_req = tau_res_u / tau_perm

# Weld leg size 's' is related to throat 't' by t = s * cos(45°)
s_req = t_req / math.cos(math.radians(45))

print("--- Step 5: Required Weld Size 's' ---")
print("The final calculation is based on the equation: s = τ_res_u / (τ_perm * cos(45°))")
print(f"s = {tau_res_u:.2f} N/mm / ({tau_perm} N/mm² * {math.cos(math.radians(45)):.3f})")
print(f"\nRequired throat thickness, t = {tau_res_u:.2f} / {tau_perm} = {t_req:.3f} mm")
print(f"Required weld leg size, s = {t_req:.3f} / {math.cos(math.radians(45)):.3f} = {s_req:.3f} mm")
print(f"\nThe determined size of the weld (leg) is {s_req:.2f} mm.")
