import math

# --- Given Data ---
P = 90000.0  # Force in N
l1 = 50.0   # Length of top weld W1 in mm
l2 = 50.0   # Length of bottom weld W2 in mm
l3 = 100.0  # Length of vertical weld W3 in mm
d1 = 50.0   # Horizontal length of the C-channel base
d2 = 150.0  # Horizontal distance from the base to the force
tau_perm = 250.0  # Permissible shear stress in N/mm^2

print("Step 1: Calculate the centroid of the weld group")
# We set the origin (0,0) at the bottom end of the vertical weld W3.
# Centroids of individual welds:
# W1: (25, 100), W2: (25, 0), W3: (0, 50)
L = l1 + l2 + l3
x_bar = (l1 * 25.0 + l2 * 25.0 + l3 * 0.0) / L
y_bar = (l1 * 100.0 + l2 * 0.0 + l3 * 50.0) / L
print(f"Total weld length, L = {L:.2f} mm")
print(f"Centroid coordinates (x_bar, y_bar) = ({x_bar:.2f} mm, {y_bar:.2f} mm)\n")


print("Step 2: Calculate the primary and secondary loads")
# Primary shear force is P
# Secondary load is the moment (torsion) M
# Eccentricity is the horizontal distance from the weld centroid to the line of action of the force
force_x_location = d1 + d2
e = force_x_location - x_bar
M = P * e
print(f"Eccentricity, e = {force_x_location:.2f} mm - {x_bar:.2f} mm = {e:.2f} mm")
print(f"Moment, M = P * e = {P:.0f} N * {e:.2f} mm = {M:.2f} N-mm\n")

print("Step 3: Calculate the polar moment of inertia of the weld group (per unit throat thickness)")
# We calculate the polar moment of inertia J_u = J / t, where t is throat thickness.
# We use the parallel axis theorem: I = I_c + A*d^2, treating welds as lines.
# Moment of inertia about the horizontal axis through the centroid (I_ux)
I_ux_1 = l1 * (100.0 - y_bar)**2
I_ux_2 = l2 * (0.0 - y_bar)**2
I_ux_3 = l3**3 / 12.0
I_ux = I_ux_1 + I_ux_2 + I_ux_3

# Moment of inertia about the vertical axis through the centroid (I_uy)
I_uy_1 = (l1**3 / 12.0) + l1 * (25.0 - x_bar)**2
I_uy_2 = (l2**3 / 12.0) + l2 * (25.0 - x_bar)**2
I_uy_3 = l3 * (0.0 - x_bar)**2
I_uy = I_uy_1 + I_uy_2 + I_uy_3

J_u = I_ux + I_uy
print(f"I_ux (unit) = {I_ux:.2f} mm^3")
print(f"I_uy (unit) = {I_uy:.2f} mm^3")
print(f"Polar Moment of Inertia, J_u = I_ux + I_uy = {J_u:.2f} mm^3\n")


print("Step 4: Calculate stress components at the critical point")
# The critical points are the corners farthest from the centroid, A(50,100) and B(50,0).
# Let's analyze point A(50, 100).
x_crit, y_crit = 50.0, 100.0
print(f"Analyzing critical point at ({x_crit:.0f}, {y_crit:.0f})")

# Primary shear stress (acts downwards, in -y direction)
# We calculate tau_prime_t = tau' * t = P / L
tau_prime_t = P / L

# Secondary shear stress components from torsion M (clockwise moment)
# tau_s_t = tau'' * t
r_x = x_crit - x_bar
r_y = y_crit - y_bar
# For a clockwise moment, τ''_x = (M/J)*r_y and τ''_y = (M/J)*(-r_x)
tau_secondary_x_t = (M / J_u) * r_y
tau_secondary_y_t = (M / J_u) * (-r_x)

# Total shear stress components
tau_total_x_t = tau_secondary_x_t
tau_total_y_t = -tau_prime_t + tau_secondary_y_t

print(f"Primary shear stress * t (τ' * t) = {P:.0f} N / {L:.0f} mm = {tau_prime_t:.2f} N/mm (acting downward)")
print(f"Secondary shear (x-comp) * t (τ''_x * t) = {tau_secondary_x_t:.2f} N/mm")
print(f"Secondary shear (y-comp) * t (τ''_y * t) = {tau_secondary_y_t:.2f} N/mm")
print(f"Total shear (x-comp) * t (τ_x * t) = {tau_total_x_t:.2f} N/mm")
print(f"Total shear (y-comp) * t (τ_y * t) = {-tau_prime_t:.2f} + {tau_secondary_y_t:.2f} = {tau_total_y_t:.2f} N/mm\n")

print("Step 5: Calculate resultant stress and required weld size")
# Resultant shear stress * t
tau_max_t = math.sqrt(tau_total_x_t**2 + tau_total_y_t**2)
print(f"Resultant shear stress * t (τ_max * t) = sqrt({tau_total_x_t:.2f}^2 + {tau_total_y_t:.2f}^2) = {tau_max_t:.2f} N/mm")

# Solve for required throat thickness t
# τ_max = τ_perm  =>  (τ_max * t) / t = τ_perm  =>  t = (τ_max * t) / τ_perm
t_req = tau_max_t / tau_perm
print(f"Required throat thickness, t = (τ_max * t) / τ_perm = {tau_max_t:.2f} N/mm / {tau_perm:.2f} N/mm^2 = {t_req:.2f} mm")

# Solve for required leg size s
# t = s * cos(45°) => s = t / cos(45°)
s_req = t_req / math.cos(math.radians(45))
print(f"Required weld leg size, s = t / cos(45°) = {t_req:.2f} mm / {math.cos(math.radians(45)):.3f} = {s_req:.2f} mm\n")

print("Final Answer:")
print(f"The required size of the weld (leg) is {s_req:.2f} mm.")
print("The final equation for the weld leg size 's' is:")
print(f"s = (sqrt({tau_total_x_t:.2f}^2 + {tau_total_y_t:.2f}^2) / {tau_perm:.0f}) / {math.cos(math.radians(45)):.3f}")
<<<17.13>>>