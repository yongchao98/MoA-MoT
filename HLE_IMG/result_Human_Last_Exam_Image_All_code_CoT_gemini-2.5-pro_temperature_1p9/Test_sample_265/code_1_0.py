import math

# --- 1. Define Constants & Given Data ---
P = 90 * 1000  # Force in N
tau_permissible = 250  # Permissible shear stress in N/mm^2
l1 = 50.0  # Length of top weld W1 in mm
l2 = 50.0  # Length of bottom weld W2 in mm
l3 = 100.0 # Length of vertical weld W3 in mm
force_arm = 150.0 # Horizontal distance from weld end to force in mm

print("This script calculates the required leg size for the welded connection.")
print("-" * 60)

# --- 2. Weld Group Properties ---
print("Step 1: Calculating Weld Group Properties")

# Coordinate system: origin at the geometric center of the vertical weld W3
# W3: runs from (0, -50) to (0, 50)
# W1: runs from (0, 50) to (50, 50)
# W2: runs from (0, -50) to (50, -50)

# Total weld length (proportional to area)
L = l1 + l2 + l3
print(f"Total weld length, L = {l1} + {l2} + {l3} = {L:.1f} mm")

# Centroid of the weld group (x_bar, y_bar)
# Sum of (Area_i * x_i) / Total Area
x_bar = (l1 * (l1 / 2) + l2 * (l2 / 2) + l3 * 0) / L
# By symmetry, y_bar is 0
y_bar = (l1 * (l3 / 2) + l2 * (-l3 / 2) + l3 * 0) / L
print(f"Centroid G = ({x_bar:.2f} mm, {y_bar:.2f} mm)")

# Eccentricity and Torque
# The force is applied at a horizontal position of l1 + force_arm from the vertical plate
x_force = l1 + force_arm
e = x_force - x_bar
T = P * e
print(f"Eccentricity, e = ({l1:.1f} + {force_arm:.1f}) - {x_bar:.2f} = {e:.2f} mm")
print(f"Torque, T = P * e = {P} N * {e:.2f} mm = {T:.2e} N-mm")

# Moment of Inertia of weld group (per unit throat thickness)
# Using parallel axis theorem: I = I_c + A*d^2 (where A is length L for unit thickness)
# For Weld 1 (W1) at y = 50mm
I_u_xx_1 = (0) + l1 * (l3/2 - y_bar)**2
I_u_yy_1 = (l1**3 / 12) + l1 * (l1/2 - x_bar)**2

# For Weld 2 (W2) at y = -50mm
I_u_xx_2 = (0) + l2 * (-l3/2 - y_bar)**2
I_u_yy_2 = (l2**3 / 12) + l2 * (l2/2 - x_bar)**2

# For Weld 3 (W3) at x = 0mm
I_u_xx_3 = (l3**3 / 12) + l3 * (0 - y_bar)**2
I_u_yy_3 = (0) + l3 * (0 - x_bar)**2

# Total Unit Moments of Inertia
I_u_xx = I_u_xx_1 + I_u_xx_2 + I_u_xx_3
I_u_yy = I_u_yy_1 + I_u_yy_2 + I_u_yy_3
# Unit Polar Moment of Inertia
J_u = I_u_xx + I_u_yy
print(f"Unit Polar Moment of Inertia, J_u = {I_u_xx:.2f} + {I_u_yy:.2f} = {J_u:.2f} mm^3\n")


# --- 3. Stress Calculation ---
print("Step 2: Calculating Stresses at Critical Point")

# Primary (Direct) Shear Stress (per unit throat thickness, denoted by prime ')
tau_p_prime = P / L
print(f"Primary shear stress (per unit throat t), τ_p' = P / L = {P} / {L:.1f} = {tau_p_prime:.2f} N/mm")

# Secondary (Torsional) Shear Stress at critical point A (top-right corner)
# Coordinates of A relative to origin: (50, 50)
# Coordinates of A relative to centroid G(12.5, 0):
r_x = l1 - x_bar
r_y = l3 / 2 - y_bar

# The stress due to a clockwise torque T has components: τ_s_x' = T*r_y/J_u and τ_s_y' = -T*r_x/J_u
tau_sx_prime = (T * r_y) / J_u
tau_sy_prime = -(T * r_x) / J_u

print("\nAnalysis at critical point A (top-right corner):")
print(f"Horizontal Torsional Stress, τ_sx' = (T * r_y) / J_u = ({T:.2e} * {r_y:.2f}) / {J_u:.2f} = {tau_sx_prime:.2f} N/mm")
print(f"Vertical Torsional Stress, τ_sy' = (-T * r_x) / J_u = (-{T:.2e} * {r_x:.2f}) / {J_u:.2f} = {tau_sy_prime:.2f} N/mm")


# Combine stresses at Point A
tau_res_x_prime = tau_sx_prime
# Primary stress (down) and torsional stress (down) add up in the vertical direction
tau_res_y_prime = tau_sy_prime - tau_p_prime

# Resultant shear stress (per unit throat)
tau_r_prime = math.sqrt(tau_res_x_prime**2 + tau_res_y_prime**2)

print("\nResultant Stress Calculation (values are per unit throat 't'):")
print(f"τ_resultant' = sqrt( (τ_sx')^2 + (τ_sy' + τ_p')^2 )")
print(f"τ_resultant' = sqrt( ({tau_sx_prime:.2f})^2 + ({tau_sy_prime:.2f} + ({-tau_p_prime:.2f}))^2 )")
print(f"τ_resultant' = sqrt( ({tau_res_x_prime:.2f})^2 + ({tau_res_y_prime:.2f})^2 ) = {tau_r_prime:.2f} N/mm\n")

# --- 4. Determine Required Weld Size ---
print("Step 3: Determining Required Weld Size")
# The actual max stress is τ_max = τ_r_prime / t
# Set τ_max = τ_permissible to find the required throat thickness t
t_req = tau_r_prime / tau_permissible
print(f"Required throat thickness, t = τ_resultant' / τ_permissible = {tau_r_prime:.2f} N/mm / {tau_permissible} N/mm^2 = {t_req:.4f} mm")

# Weld leg size 's' is related to throat 't' by t = 0.707 * s
s_req = t_req / 0.707
print(f"For a standard fillet weld, leg size s = t / 0.707")
print(f"Required leg size, s = {t_req:.4f} / 0.707 = {s_req:.2f} mm")
print("-" * 60)
print(f"The required size of the welds (leg) should be at least {s_req:.2f} mm.")