import math

# --- 1. Given Data ---
P = 90 * 1000  # Force in N
tau_perm = 250  # Permissible shear stress in N/mm^2 (MPa)

# Weld dimensions in mm
l1 = 50.0  # Length of top weld W1
l2 = 50.0  # Length of bottom weld W2
l3 = 100.0 # Length of vertical weld W3
e_prime = 150.0 # Distance from weld end to force

# --- 2. Weld Group Properties ---

# Total length of the weld
L_total = l1 + l2 + l3

# Centroid calculation (origin at the intersection of W2 and W3)
# x_bar = (A1*x1 + A2*x2 + ...) / (A1 + A2 + ...)
# Here, we use length L instead of area A since throat thickness is constant
x_bar = (l1 * (l1 / 2) + l2 * (l2 / 2) + l3 * 0) / L_total
y_bar = (l1 * l3 + l2 * 0 + l3 * (l3 / 2)) / L_total

# Eccentricity 'e' and Torque 'T'
# The force is applied at a distance of (l1 + e_prime) from the vertical weld line.
x_force = l1 + e_prime
e = x_force - x_bar
T = P * e

# Unit Moment of Inertia calculation about the centroid
# Using parallel axis theorem: I = I_c + A*d^2 (using L for A)
# I_xx (about centroidal x-axis at y_bar)
I_xx_1 = (l1 * (l3 - y_bar)**2)
I_xx_2 = (l2 * (0 - y_bar)**2)
I_xx_3 = (l3**3 / 12) + l3 * ((l3/2) - y_bar)**2
I_xx_u = I_xx_1 + I_xx_2 + I_xx_3

# I_yy (about centroidal y-axis at x_bar)
I_yy_1 = (l1**3 / 12) + l1 * ((l1/2) - x_bar)**2
I_yy_2 = (l2**3 / 12) + l2 * ((l2/2) - x_bar)**2
I_yy_3 = l3 * (0 - x_bar)**2
I_yy_u = I_yy_1 + I_yy_2 + I_yy_3

# Unit Polar Moment of Inertia
J_u = I_xx_u + I_yy_u

# --- 3. Stress Analysis at Critical Point ---
# The maximum stress occurs at the corners farthest from the centroid.
# Let's analyze the top-right corner of the weld (end of W1).
# Coordinates of this point relative to the centroid:
x_r = l1 - x_bar
y_r = l3 - y_bar

# Primary shear stress components (acting downwards)
# tau_p = P / (L_total * t)
tau_p_y_per_t = -P / L_total  # Stress per unit throat thickness 't'
tau_p_x_per_t = 0

# Secondary (torsional) shear stress components
# tau_s = (T*r)/J. J = J_u * t.
# Components: tau_sx = -T*y / J, tau_sy = T*x / J
tau_s_x_per_t = (-T * y_r) / J_u
tau_s_y_per_t = (T * x_r) / J_u

# Total shear stress components
tau_x_total_per_t = tau_p_x_per_t + tau_s_x_per_t
tau_y_total_per_t = tau_p_y_per_t + tau_s_y_per_t

# --- 4. Determine Weld Size ---
# The resultant stress magnitude is tau_res = sqrt(tau_x^2 + tau_y^2)
# tau_perm = (1/t) * sqrt(tau_x_total_per_t^2 + tau_y_total_per_t^2)
# Solve for throat thickness 't'
t = (1 / tau_perm) * math.sqrt(tau_x_total_per_t**2 + tau_y_total_per_t**2)

# Convert throat thickness 't' to leg size 's' (t = s * cos(45) = s * 0.707)
s = t / (math.cos(math.radians(45)))

# --- 5. Print Results ---
print("--- Analysis of Welded Connection ---")
print(f"Force, P = {P/1000} kN")
print(f"Permissible Shear Stress, τ_perm = {tau_perm} N/mm²\n")

print("1. Weld Group Properties:")
print(f"  - Total Weld Length, L = {L_total:.2f} mm")
print(f"  - Centroid (x_bar, y_bar) = ({x_bar:.2f}, {y_bar:.2f}) mm")
print(f"  - Eccentricity, e = {e:.2f} mm")
print(f"  - Torque, T = P * e = {T:.2f} N-mm")
print(f"  - Unit Polar Moment of Inertia, J_u = {J_u:.2f} mm³\n")

print("2. Stresses at Critical Point (Top-Right Corner):")
print(f"  - Primary Shear (Vertical), τ_p' = P / L = {abs(tau_p_y_per_t):.2f} / t  N/mm²")
print(f"  - Secondary Shear (Horizontal), τ_sx' = -(T*y_r)/J_u = {tau_s_x_per_t:.2f} / t N/mm²")
print(f"  - Secondary Shear (Vertical), τ_sy' = (T*x_r)/J_u = {tau_s_y_per_t:.2f} / t N/mm²\n")

print("3. Resultant Stress Calculation:")
print("   τ_res = √((τ_sx')² + (τ_p' + τ_sy')²)")
final_eq = (f"   {tau_perm} = (1/t) * √(({tau_s_x_per_t:.2f})² "
            f"+ ({tau_p_y_per_t:.2f} + {tau_s_y_per_t:.2f})²)")
print(final_eq)

print("\n4. Required Weld Size:")
print(f"Solving the equation gives the required throat thickness, t = {t:.4f} mm")
print("The leg size 's' is calculated from t = s * cos(45°).")
print(f"Required Leg Size, s = t / 0.707 = {s:.4f} mm")
print("\nFinal Answer (rounded to one decimal place):")
print(f"The required size of the welds (leg) is {s:.1f} mm.")
<<<14.1>>>