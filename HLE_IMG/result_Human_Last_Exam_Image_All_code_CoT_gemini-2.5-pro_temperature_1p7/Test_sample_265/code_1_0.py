import math

# --- 1. Given Data ---
P = 90000.0  # Force in N (90 kN)
tau_perm = 250.0  # Permissible shear stress in N/mm^2
# Dimensions from the image
weld_length_horiz = 50.0  # Length of horizontal welds (W1, W2) in mm
weld_length_vert = 100.0  # Length of the vertical weld (W3) in mm
ecc_dist_from_vert_weld = 150.0 # Distance from W3 to the force line of action

print("--- Step 1: Find the Centroid of the Weld Group ---")
# Weld segments as (length, x_centroid, y_centroid) with origin at bottom-left corner of W3
w1 = (weld_length_horiz, weld_length_horiz / 2, weld_length_vert)
w2 = (weld_length_horiz, weld_length_horiz / 2, 0)
w3 = (weld_length_vert, 0, weld_length_vert / 2)
welds = [w1, w2, w3]

total_length = sum(w[0] for w in welds)
sum_lx = sum(w[0] * w[1] for w in welds)
sum_ly = sum(w[0] * w[2] for w in welds)

x_bar = sum_lx / total_length
y_bar = sum_ly / total_length

print(f"Total weld length (L) = {total_length:.2f} mm")
print(f"Centroid coordinates (x_bar, y_bar) = ({x_bar:.2f} mm, {y_bar:.2f} mm)")

print("\n--- Step 2: Calculate Unit Polar Moment of Inertia (J_u) ---")
# J_u is the polar moment of inertia for a unit throat thickness (t=1 mm)

# I_xx for horizontal welds (W1, W2)
I_xx_1 = (w1[0] * (w1[2] - y_bar)**2)
I_xx_2 = (w2[0] * (w2[2] - y_bar)**2)
# I_xx for vertical weld (W3)
I_xx_3 = (w3[0]**3 / 12) + w3[0] * (w3[2] - y_bar)**2
I_xx_u = I_xx_1 + I_xx_2 + I_xx_3

# I_yy for horizontal welds (W1, W2)
I_yy_1 = (w1[0]**3 / 12) + w1[0] * (w1[1] - x_bar)**2
I_yy_2 = (w2[0]**3 / 12) + w2[0] * (w2[1] - x_bar)**2
# I_yy for vertical weld (W3)
I_yy_3 = w3[0] * (w3[1] - x_bar)**2
I_yy_u = I_yy_1 + I_yy_2 + I_yy_3

J_u = I_xx_u + I_yy_u
print(f"Unit Moment of Inertia about x-axis (I_xx_u) = {I_xx_u:.2f} mm^3")
print(f"Unit Moment of Inertia about y-axis (I_yy_u) = {I_yy_u:.2f} mm^3")
print(f"Unit Polar Moment of Inertia (J_u) = {J_u:.2f} mm^3")

print("\n--- Step 3: Calculate Stresses at the Critical Point ---")
# The critical points are the corners farthest from the centroid.
# Let's analyze the top-right corner, point A at (50, 100).
crit_point_x = weld_length_horiz
crit_point_y = weld_length_vert

# Eccentricity and Moment
force_x_location = weld_length_horiz + ecc_dist_from_vert_weld
e = force_x_location - x_bar
M = P * e
print(f"Eccentricity (e) = {e:.2f} mm")
print(f"Bending Moment (M) = {M:.2f} N-mm")

# Primary Shear Stress (per unit throat 't')
# Acts downwards (-y direction)
tau_prime_per_t = P / total_length
print(f"Primary Shear Stress (τ') = {tau_prime_per_t:.2f} / t N/mm^2")

# Secondary Shear Stress at the critical point (per unit throat 't')
r_x = crit_point_x - x_bar
r_y = crit_point_y - y_bar
r = math.sqrt(r_x**2 + r_y**2)

# Secondary stress components due to torsion
# Moment is clockwise, resisting stress is counter-clockwise
tau_second_x_per_t = - (M * r_y) / J_u
tau_second_y_per_t = (M * r_x) / J_u
tau_second_mag_per_t = math.sqrt(tau_second_x_per_t**2 + tau_second_y_per_t**2)

print(f"Analyzing critical point A at ({crit_point_x}, {crit_point_y}) which is at radius r = {r:.2f} mm")
print(f"Secondary Shear Stress Components (τ''):")
print(f"  τ''_x = {tau_second_x_per_t:.2f} / t N/mm^2")
print(f"  τ''_y = {tau_second_y_per_t:.2f} / t N/mm^2")

print("\n--- Step 4: Find Maximum Resultant Stress ---")
# Total stress components (vector sum)
total_tau_x_per_t = 0 + tau_second_x_per_t
total_tau_y_per_t = -tau_prime_per_t + tau_second_y_per_t

resultant_tau_per_t = math.sqrt(total_tau_x_per_t**2 + total_tau_y_per_t**2)

print(f"Total x-component τ_x = {total_tau_x_per_t:.2f} / t")
print(f"Total y-component τ_y = -{tau_prime_per_t:.2f}/t + {tau_second_y_per_t:.2f}/t = {total_tau_y_per_t:.2f} / t")
print(f"Resultant Shear Stress (τ_max) = {resultant_tau_per_t:.2f} / t")


print("\n--- Step 5: Determine the Weld Size (s) ---")
# Equate resultant stress with permissible stress: τ_max = τ_perm
# resultant_tau_per_t / t = tau_perm
t = resultant_tau_per_t / tau_perm
# For a fillet weld, t = s * cos(45deg) = s * 0.707
s = t / (math.cos(math.radians(45)))

print("\nThe final equation is τ_max = τ_perm")
print(f"{resultant_tau_per_t:.2f} / t = {tau_perm:.2f}")
print(f"Required throat thickness (t) = {resultant_tau_per_t:.2f} / {tau_perm:.2f} = {t:.4f} mm")
print("\nConvert throat thickness (t) to leg size (s) using s = t / 0.707")
print(f"s = {t:.4f} / {math.cos(math.radians(45)):.4f}")
print(f"Required weld leg size (s) = {s:.2f} mm")
print("\nFinal Answer:")
print(round(s, 2))