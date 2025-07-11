import math

# --- Given Data ---
P = 90 * 1000  # Force in N
tau_perm = 250  # Permissible shear stress in N/mm^2
l1 = 50.0       # Length of weld W1 in mm
l2 = 50.0       # Length of weld W2 in mm
l3 = 100.0      # Length of weld W3 in mm
ecc_dist = 150.0 # Eccentric distance from weld W3 in mm

print("--- Step 1: Find the Centroid of the Weld Group ---")
# We set the origin (0,0) at the bottom end of the vertical weld W3.
# Centroid of W1: (25, 100)
# Centroid of W2: (25, 0)
# Centroid of W3: (0, 50)

L_total = l1 + l2 + l3
x_bar = (l1 * 25 + l2 * 25 + l3 * 0) / L_total
y_bar = (l1 * 100 + l2 * 0 + l3 * 50) / L_total

print(f"Total length of weld, L = {l1} + {l2} + {l3} = {L_total} mm")
print(f"x_bar = ({l1}*25 + {l2}*25 + {l3}*0) / {L_total} = {x_bar:.2f} mm")
print(f"y_bar = ({l1}*100 + {l2}*0 + {l3}*50) / {L_total} = {y_bar:.2f} mm")
print(f"The centroid G is at ({x_bar:.2f}, {y_bar:.2f}) mm.\n")

print("--- Step 2: Calculate Moments of Inertia (per unit throat thickness) ---")
# Using Parallel Axis Theorem: I = I_c + A*d^2
# For a line, I_c about its strong axis is l^3/12 and negligible about its weak axis.
# All calculations are per unit throat thickness 't', so Area 'A' is replaced by length 'l'.

# Ixx about the centroidal axis (y = 50)
dy1 = 100 - y_bar
dy2 = 0 - y_bar
dy3 = 50 - y_bar
Ixx_unit = (l1 * dy1**2) + (l2 * dy2**2) + (l3**3 / 12 + l3 * dy3**2)
print(f"I_xx = (l1*dy1^2) + (l2*dy2^2) + (l3^3/12) = ({l1}*{dy1**2:.0f}) + ({l2}*{dy2**2:.0f}) + ({l3**3 / 12:.0f}) = {Ixx_unit:.2f} mm^3")

# Iyy about the centroidal axis (x = 12.5)
dx1 = 25 - x_bar
dx2 = 25 - x_bar
dx3 = 0 - x_bar
Iyy_unit = (l1**3 / 12 + l1 * dx1**2) + (l2**3 / 12 + l2 * dx2**2) + (l3 * dx3**2)
print(f"I_yy = (l1^3/12 + l1*dx1^2) + (l2^3/12 + l2*dx2^2) + (l3*dx3^2) = {Iyy_unit:.2f} mm^3")

# Polar Moment of Inertia
J_unit = Ixx_unit + Iyy_unit
print(f"Polar Moment of Inertia, J = I_xx + I_yy = {Ixx_unit:.2f} + {Iyy_unit:.2f} = {J_unit:.2f} mm^3\n")


print("--- Step 3: Calculate Stresses at the Critical Point ---")
# The critical points are the corners farthest from the centroid: A(50, 100) and B(50, 0).
# Let's analyze point A(50, 100).
x_A = 50
y_A = 100

# Eccentricity and Torque
e = (l1 + ecc_dist) - x_bar
T = P * e
print(f"Eccentricity, e = ({l1} + {ecc_dist}) - {x_bar:.2f} = {e:.2f} mm")
print(f"Torque, T = P * e = {P} N * {e:.2f} mm = {T:.2f} N-mm\n")

# Primary Shear Stress (τ') - acts vertically downwards
# Note: We calculate stress multiplied by throat thickness 't'
tau_prime_y_unit = -P / L_total
print(f"Primary shear stress (τ') component (multiplied by t):")
print(f"τ'_y * t = -P / L = -{P} / {L_total} = {tau_prime_y_unit:.2f} N/mm")

# Secondary Shear Stress (τ'') - due to torque
# Position of point A relative to centroid
r_x = x_A - x_bar
r_y = y_A - y_bar
# Components of secondary shear stress (multiplied by t)
tau_double_prime_x_unit = (T * r_y) / J_unit
tau_double_prime_y_unit = -(T * r_x) / J_unit
print(f"\nSecondary shear stress (τ'') components at point A({x_A},{y_A}) (multiplied by t):")
print(f"τ''_x * t = (T * r_y) / J = ({T:.0f} * {r_y:.2f}) / {J_unit:.2f} = {tau_double_prime_x_unit:.2f} N/mm")
print(f"τ''_y * t = -(T * r_x) / J = -({T:.0f} * {r_x:.2f}) / {J_unit:.2f} = {tau_double_prime_y_unit:.2f} N/mm")

# Resultant Shear Stress (τ_res)
tau_res_x_unit = tau_double_prime_x_unit
tau_res_y_unit = tau_prime_y_unit + tau_double_prime_y_unit
print(f"\nResultant stress components (multiplied by t):")
print(f"τ_x * t = {tau_res_x_unit:.2f} N/mm")
print(f"τ_y * t = {tau_prime_y_unit:.2f} + {tau_double_prime_y_unit:.2f} = {tau_res_y_unit:.2f} N/mm")

# Magnitude of resultant stress
tau_res_unit = math.sqrt(tau_res_x_unit**2 + tau_res_y_unit**2)
print(f"\nResultant stress magnitude, τ_max = (sqrt( ({tau_res_x_unit:.2f})^2 + ({tau_res_y_unit:.2f})^2 )) / t")
print(f"τ_max = {tau_res_unit:.2f} / t\n")

print("--- Step 4: Determine the Required Weld Size ---")
# Set τ_max = τ_perm to solve for the throat thickness 't'
# τ_perm = tau_res_unit / t  =>  t = tau_res_unit / τ_perm
t = tau_res_unit / tau_perm
print(f"Required throat thickness, t = τ_max_unit / τ_perm = {tau_res_unit:.2f} / {tau_perm} = {t:.2f} mm")

# Convert throat thickness 't' to leg size 's' for a fillet weld
s = t / 0.707
print(f"Required weld leg size, s = t / 0.707 = {t:.2f} / 0.707 = {s:.2f} mm")

print("\nFinal Answer:")
print(f"The required size of the welds (leg) is {s:.2f} mm.")
print(f"<<<{s:.2f}>>>")