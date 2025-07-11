import math

# --- 1. Define Given Data ---
# All dimensions are in mm, forces in N, stress in N/mm^2
P = 90000.0  # Force
tau_perm = 250.0  # Permissible shear stress
l1 = 50.0  # Length of top horizontal weld (W1)
l2 = 50.0  # Length of bottom horizontal weld (W2)
l3 = 100.0 # Length of vertical weld (W3)
e_arm = 150.0 # Eccentricity arm from the end of the weld

# --- 2. Calculate Weld Group Properties ---
print("--- Step 1: Weld Group Properties ---")
# Total length of the weld
L_total = l1 + l2 + l3

# Centroid Calculation (origin at the bottom of the vertical weld W3)
# Centroid of W1: (25, 100), Centroid of W2: (25, 0), Centroid of W3: (0, 50)
x_bar = (l1 * 25.0 + l2 * 25.0 + l3 * 0.0) / L_total
y_bar = (l1 * 100.0 + l2 * 0.0 + l3 * 50.0) / L_total

print(f"Total weld length, L = {l1} + {l2} + {l3} = {L_total:.2f} mm")
print(f"The centroid of the weld group (x_bar, y_bar) is ({x_bar:.2f} mm, {y_bar:.2f} mm).")
print("-" * 40)

# --- 3. Calculate Primary Shear Stress ---
print("--- Step 2: Primary Shear Stress (τ') ---")
# This is the stress per unit throat thickness (t)
tau_prime_per_t = P / L_total
print("Primary shear stress is due to the direct force P.")
print(f"τ' / t = P / L_total = {P:.2f} / {L_total:.2f} = {tau_prime_per_t:.2f} N/mm^2 per mm of throat thickness.")
print("-" * 40)


# --- 4. Calculate Secondary Shear Stress ---
print("--- Step 3: Secondary Shear Stress (τ'') ---")
# Moment of Inertia per unit throat thickness (J_u)
# I_xx for the group
I_xx_u = (l1 * (100.0 - y_bar)**2) + (l2 * (0.0 - y_bar)**2) + (l3**3 / 12)
# I_yy for the group
I_yy_u = (l1**3 / 12 + l1 * (25.0 - x_bar)**2) + \
         (l2**3 / 12 + l2 * (25.0 - x_bar)**2) + \
         (l3 * (0.0 - x_bar)**2)
# Polar moment of inertia per unit throat thickness
J_u = I_xx_u + I_yy_u

# Eccentricity and Torque
e = (l1 + e_arm) - x_bar
T = P * e
print("Secondary shear stress is due to the torque T = P * e.")
print(f"Polar Moment of Inertia, J_u = I_xx_u + I_yy_u = {I_xx_u:.2f} + {I_yy_u:.2f} = {J_u:.2f} mm^3")
print(f"Eccentricity, e = ({l1:.2f} + {e_arm:.2f}) - {x_bar:.2f} = {e:.2f} mm")
print(f"Torque, T = {P:.2f} N * {e:.2f} mm = {T:.2f} N-mm")
print("-" * 40)

# --- 5. Combine Stresses at Critical Point ---
print("--- Step 4: Resultant Stress at Critical Point ---")
# Critical points are the top-right (50, 100) and bottom-right (50, 0) corners.
# Let's analyze the bottom-right corner: B(50, 0)
x_rel = 50.0 - x_bar
y_rel = 0.0 - y_bar
print(f"Analyzing critical point B({50.0:.2f}, {0.0:.2f}), which is furthest from the centroid.")
print(f"Relative position of B to centroid: x_rel = {x_rel:.2f} mm, y_rel = {y_rel:.2f} mm")

# Stress components per unit throat thickness
# Primary (acts downward)
tau_prime_y_per_t = -tau_prime_per_t
# Secondary (from clockwise torque)
tau_double_prime_x_per_t = (T * y_rel) / J_u
tau_double_prime_y_per_t = (-T * x_rel) / J_u

# Total stress components
tau_total_x_per_t = tau_double_prime_x_per_t
tau_total_y_per_t = tau_prime_y_per_t + tau_double_prime_y_per_t

print("\nCombining stress components:")
print(f"τ_x/t = (T * y_rel) / J_u = ({T:.2f} * {y_rel:.2f}) / {J_u:.2f} = {tau_total_x_per_t:.2f}")
print(f"τ_y/t = (P/L) + (-T * x_rel) / J_u = {tau_prime_y_per_t:.2f} + {tau_double_prime_y_per_t:.2f} = {tau_total_y_per_t:.2f}")

# Resultant stress factor K where τ_res = K / t
K = math.sqrt(tau_total_x_per_t**2 + tau_total_y_per_t**2)

print("\nResultant stress equation:")
print(f"τ_res = (1/t) * sqrt( (τ_x)^2 + (τ_y)^2 )")
print(f"{tau_perm:.2f} = (1/t) * sqrt( ({tau_total_x_per_t:.2f})^2 + ({tau_total_y_per_t:.2f})^2 )")
print(f"{tau_perm:.2f} = (1/t) * {K:.2f}")
print("-" * 40)

# --- 6. Determine Weld Size ---
print("--- Step 5: Calculate Required Weld Size (s) ---")
# Required throat thickness t
t_req = K / tau_perm
# Required weld leg size s
s_req = t_req / 0.707

print(f"Required throat thickness, t = {K:.2f} / {tau_perm:.2f} = {t_req:.2f} mm")
print(f"Required weld size (leg), s = t / 0.707 = {t_req:.2f} / 0.707 = {s_req:.2f} mm")
print("-" * 40)

print(f"\nThe final calculated size of the welds (leg) is {s_req:.2f} mm.")
<<<17.13>>>