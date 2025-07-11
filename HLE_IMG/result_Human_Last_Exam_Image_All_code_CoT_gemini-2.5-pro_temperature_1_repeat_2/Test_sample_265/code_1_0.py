import math

# --- 1. Given Data ---
P = 90 * 1000  # Force in N
tau_permissible = 250  # Permissible shear stress in N/mm^2
l1 = 50  # Length of weld W1 in mm
l2 = 50  # Length of weld W2 in mm
l3 = 100 # Length of weld W3 in mm
e_dist = 150 # Eccentric distance from weld W1/W2 to force in mm

print("--- Step 1: Define Input Parameters ---")
print(f"Force, P = {P/1000} kN = {P} N")
print(f"Permissible Shear Stress, τ_perm = {tau_permissible} N/mm^2")
print(f"Weld W1 length, l1 = {l1} mm")
print(f"Weld W2 length, l2 = {l2} mm")
print(f"Weld W3 length, l3 = {l3} mm")
print(f"Eccentricity distance from vertical plate, e_dist = {e_dist} mm\n")

# --- 2. Centroid Calculation ---
# Let the origin (0,0) be at the bottom of the vertical weld W3.
# Centroid of W1: (x1, y1) = (25, 100)
# Centroid of W2: (x2, y2) = (25, 0)
# Centroid of W3: (x3, y3) = (0, 50)
x1, y1 = l1 / 2, l3
x2, y2 = l2 / 2, 0
x3, y3 = 0, l3 / 2

L = l1 + l2 + l3 # Total length of the weld
x_bar = (l1 * x1 + l2 * x2 + l3 * x3) / L
y_bar = (l1 * y1 + l2 * y2 + l3 * y3) / L

print("--- Step 2: Calculate the Centroid of the Weld Group G(x_bar, y_bar) ---")
print("Using origin at the bottom-left corner of the weld pattern.")
print(f"Total weld length, L = l1 + l2 + l3 = {l1} + {l2} + {l3} = {L} mm")
print(f"x_bar = (l1*x1 + l2*x2 + l3*x3) / L = ({l1}*{x1:.1f} + {l2}*{x2:.1f} + {l3}*{x3:.1f}) / {L} = {x_bar:.2f} mm")
print(f"y_bar = (l1*y1 + l2*y2 + l3*y3) / L = ({l1}*{y1:.1f} + {l2}*{y2:.1f} + {l3}*{y3:.1f}) / {L} = {y_bar:.2f} mm\n")

# --- 3. Moment Calculation ---
# Eccentricity is the distance from the force to the centroid.
e = (l1 + e_dist) - x_bar
M = P * e

print("--- Step 3: Calculate the Moment (M) about the Centroid ---")
print(f"Eccentricity, e = (l1 + e_dist) - x_bar = ({l1} + {e_dist}) - {x_bar:.2f} = {e:.2f} mm")
print(f"Moment, M = P * e = {P} N * {e:.2f} mm = {M:.2f} N-mm\n")

# --- 4. Polar Moment of Inertia Calculation (per unit throat thickness) ---
# Second moment of area about x-axis (I_ux)
I_ux_1 = l1 * (y1 - y_bar)**2
I_ux_2 = l2 * (y2 - y_bar)**2
I_ux_3 = (l3**3 / 12) + l3 * (y3 - y_bar)**2
I_ux = I_ux_1 + I_ux_2 + I_ux_3

# Second moment of area about y-axis (I_uy)
I_uy_1 = (l1**3 / 12) + l1 * (x1 - x_bar)**2
I_uy_2 = (l2**3 / 12) + l2 * (x2 - x_bar)**2
I_uy_3 = l3 * (x3 - x_bar)**2
I_uy = I_uy_1 + I_uy_2 + I_uy_3

# Unit Polar Moment of Inertia
J_u = I_ux + I_uy

print("--- Step 4: Calculate the Unit Polar Moment of Inertia (J_u) ---")
print(f"I_ux = (l3^3/12 + l3*dy3^2) + (l1*dy1^2) + (l2*dy2^2)")
print(f"I_ux = ({l3**3 / 12:.2f} + {l3}*({y3-y_bar:.2f})^2) + ({l1}*({y1-y_bar:.2f})^2) + ({l2}*({y2-y_bar:.2f})^2) = {I_ux:.2f} mm^3")
print(f"I_uy = (I_uy_w1) + (I_uy_w2) + (I_uy_w3)")
print(f"I_uy = (l1^3/12 + l1*dx1^2) + (l2^3/12 + l2*dx2^2) + (l3*dx3^2) = {I_uy:.2f} mm^3")
print(f"J_u = I_ux + I_uy = {I_ux:.2f} + {I_uy:.2f} = {J_u:.2f} mm^3\n")


# --- 5. Stress Analysis at Critical Point ---
# The critical points are the corners farthest from the centroid (0,0) and (0,100).
# Let's analyze the point at (0,0) which is one of the most stressed points.
# Coordinates of critical point relative to centroid:
x_c = 0 - x_bar
y_c = 0 - y_bar

# Primary shear stress (acts in -y direction) per unit throat thickness
tau_prime_y_per_t = -P / L

# Secondary shear stress components per unit throat thickness
tau_sec_x_per_t = -M * y_c / J_u
tau_sec_y_per_t = M * x_c / J_u

# Resultant stress components per unit throat thickness
tau_res_x_per_t = tau_sec_x_per_t
tau_res_y_per_t = tau_prime_y_per_t + tau_sec_y_per_t

# Magnitude of resultant stress per unit throat thickness
tau_res_per_t = math.sqrt(tau_res_x_per_t**2 + tau_res_y_per_t**2)

print("--- Step 5: Calculate Stresses at Critical Point (0,0) ---")
print("Stresses are first calculated 'per unit throat thickness (t)'.")
print(f"Primary Shear Stress: τ'_y = P / L = {P} / {L} = {abs(tau_prime_y_per_t):.2f} / t  (downwards)")
print(f"Secondary Shear Stress (x-comp): τ''_x = -M*y / J_u = -({M:.0f})*({y_c:.2f}) / {J_u:.0f} = {tau_sec_x_per_t:.2f} / t")
print(f"Secondary Shear Stress (y-comp): τ''_y = M*x / J_u = ({M:.0f})*({x_c:.2f}) / {J_u:.0f} = {tau_sec_y_per_t:.2f} / t")
print(f"Resultant Shear (x-comp): τ_x = {tau_res_x_per_t:.2f} / t")
print(f"Resultant Shear (y-comp): τ_y = τ'_y + τ''_y = {tau_prime_y_per_t:.2f} + {tau_sec_y_per_t:.2f} = {tau_res_y_per_t:.2f} / t")
print(f"Resultant Shear Magnitude: τ_res = sqrt(τ_x^2 + τ_y^2) = sqrt(({tau_res_x_per_t:.2f}/t)^2 + ({tau_res_y_per_t:.2f}/t)^2) = {tau_res_per_t:.2f} / t\n")

# --- 6. Solve for Weld Size ---
# τ_permissible = τ_res = tau_res_per_t / t
t = tau_res_per_t / tau_permissible
s = t / 0.707

print("--- Step 6: Determine Required Weld Size (s) ---")
print("Set the maximum resultant stress equal to the permissible stress to find the throat thickness (t).")
print(f"τ_perm = τ_res => {tau_permissible} = {tau_res_per_t:.2f} / t")
print(f"Required throat thickness, t = {tau_res_per_t:.2f} / {tau_permissible} = {t:.3f} mm")
print("\nFinally, calculate the weld leg size (s) from the throat thickness.")
print(f"Weld leg size, s = t / 0.707 = {t:.3f} / 0.707 = {s:.2f} mm")

# Round the final answer to one decimal place as is conventional
s_final = round(s, 1)
print(f"\nThe required weld leg size is approximately {s_final} mm.")
print(f"<<<{s_final}>>>")
