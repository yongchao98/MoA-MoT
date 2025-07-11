import math

# --- Given Data ---
P = 90 * 1000  # Force in N (90 kN)
tau_perm = 250  # Permissible shear stress in N/mm^2
l1 = 50.0  # Length of top horizontal weld (W1) in mm
l2 = 50.0  # Length of bottom horizontal weld (W2) in mm
l3 = 100.0 # Length of vertical weld (W3) in mm
e_arm = 150.0 # Distance from the end of the beam to the force application point

# --- Step 1: Calculate Weld Properties ---
# Total length of the weld
L = l1 + l2 + l3

# Centroid Calculation (origin at bottom of W3)
# Centroids of individual welds: W1(25, 100), W2(25, 0), W3(0, 50)
x_bar = (l1 * 25 + l2 * 25 + l3 * 0) / L
y_bar = (l1 * 100 + l2 * 0 + l3 * 50) / L

print("--- Weld Group Properties ---")
print(f"Total weld length, L = {l1} + {l2} + {l3} = {L:.2f} mm")
print(f"Centroid of the weld group, G = ({x_bar:.2f} mm, {y_bar:.2f} mm)")

# --- Step 2: Calculate Moment and Primary Shear ---
# Eccentricity (distance from force to centroid)
e = (l1 + e_arm) - x_bar
# Moment (Torque)
M = P * e
# Primary shear stress (per unit throat thickness)
tau_p_per_t = P / L

print("\n--- Force and Moment Analysis ---")
print(f"Eccentricity, e = ({l1} + {e_arm}) - {x_bar:.2f} = {e:.2f} mm")
print(f"Moment, M = P * e = {P} N * {e:.2f} mm = {M:.2f} N-mm")
print(f"Primary shear force per unit throat thickness, P/A = P/(L*t) = {P}/{L}*t = {tau_p_per_t:.2f}/t N/mm^2")

# --- Step 3: Calculate Polar Moment of Inertia (J_u) ---
# Moment of inertia about x-axis (I_xx_u)
# Using parallel axis theorem: I = I_c + A*d^2 (here A is length l)
I_xx1 = (0) + l1 * (100 - y_bar)**2
I_xx2 = (0) + l2 * (0 - y_bar)**2
I_xx3 = (l3**3 / 12) + l3 * (50 - y_bar)**2
I_xx_u = I_xx1 + I_xx2 + I_xx3

# Moment of inertia about y-axis (I_yy_u)
I_yy1 = (l1**3 / 12) + l1 * (25 - x_bar)**2
I_yy2 = (l2**3 / 12) + l2 * (25 - x_bar)**2
I_yy3 = (0) + l3 * (0 - x_bar)**2
I_yy_u = I_yy1 + I_yy2 + I_yy3

# Unit polar moment of inertia
J_u = I_xx_u + I_yy_u

print("\n--- Moment of Inertia Calculation (per unit throat thickness) ---")
print(f"I_xx_u = {I_xx_u:.2f} mm^3")
print(f"I_yy_u = {I_yy_u:.2f} mm^3")
print(f"Unit Polar Moment of Inertia, J_u = I_xx_u + I_yy_u = {J_u:.2f} mm^3")

# --- Step 4: Calculate Stresses at the Critical Point ---
# The critical point is the one farthest from the centroid, which is the top-right corner (50, 100).
x_crit = 50.0
y_crit = 100.0
r_x = x_crit - x_bar
r_y = y_crit - y_bar

# Secondary shear stress components (per unit throat thickness)
# τ_s = (M*r)/J = (M*r)/(J_u*t)
tau_sx_per_t = (M * r_y) / J_u  # Horizontal component
tau_sy_per_t = (M * r_x) / J_u  # Vertical component (due to clockwise moment)

# Total stress components (per unit throat thickness)
tau_total_x_per_t = tau_sx_per_t
tau_total_y_per_t = tau_p_per_t + tau_sy_per_t

print("\n--- Stress Calculation at Critical Point ({x_crit}, {y_crit}) ---")
print(f"Primary shear stress (vertical), τ_p = {tau_p_per_t:.2f}/t")
print(f"Secondary shear stress (horizontal), τ_sx = (M*r_y)/J_u = ({M:.0f}*{r_y:.2f})/{J_u:.2f} / t = {tau_sx_per_t:.2f}/t")
print(f"Secondary shear stress (vertical), τ_sy = (M*r_x)/J_u = ({M:.0f}*{r_x:.2f})/{J_u:.2f} / t = {tau_sy_per_t:.2f}/t")

# --- Step 5: Calculate Resultant Stress and Required Weld Size ---
# Resultant stress τ_res = sqrt(τ_total_x^2 + τ_total_y^2)
# τ_res = (1/t) * sqrt(tau_total_x_per_t^2 + tau_total_y_per_t^2)
res_stress_per_t = math.sqrt(tau_total_x_per_t**2 + tau_total_y_per_t**2)

# Equate to permissible stress: τ_res = τ_perm
# res_stress_per_t / t = tau_perm  => t = res_stress_per_t / tau_perm
t = res_stress_per_t / tau_perm

# Calculate leg size s from throat thickness t (t = 0.707 * s)
s = t / 0.707

print("\n--- Final Calculation ---")
print("The resultant shear stress τ_res is the vector sum of the primary and secondary shear stresses.")
print(f"τ_res = sqrt( (τ_p + τ_sy)^2 + τ_sx^2 )")
print(f"τ_res = sqrt( ({tau_p_per_t:.2f}/t + {tau_sy_per_t:.2f}/t)^2 + ({tau_sx_per_t:.2f}/t)^2 )")
print(f"τ_res = sqrt( ({-tau_total_y_per_t:.2f}/t)^2 + ({tau_total_x_per_t:.2f}/t)^2 )")
print(f"τ_res = {res_stress_per_t:.2f} / t")
print("\nSet the resultant stress equal to the permissible stress:")
print(f"τ_res = τ_perm => {res_stress_per_t:.2f} / t = {tau_perm}")
print(f"Required throat thickness, t = {res_stress_per_t:.2f} / {tau_perm} = {t:.2f} mm")
print("\nThe leg size 's' is calculated from the throat thickness 't' (t = 0.707 * s):")
print(f"Required leg size, s = t / 0.707 = {t:.2f} / 0.707 = {s:.2f} mm")

print(f"\nThe required size of the weld (leg) is {s:.2f} mm.")
print(f"\n<<<The required size of the weld (leg) is {s:.2f} mm.>>>")
<<<17.13>>>