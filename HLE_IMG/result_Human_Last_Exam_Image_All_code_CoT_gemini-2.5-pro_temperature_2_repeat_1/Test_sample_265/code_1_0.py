import math

# --- 1. GIVEN DATA & GEOMETRY ---
P = 90000  # Force in N
tau_perm = 250  # Permissible shear stress in N/mm^2

# Dimensions from the image (in mm)
w1_length = 50.0
w2_length = 50.0
w3_length = 100.0
eccentricity_arm = 150.0

print("Step 1: Calculating the Centroid of the Weld Group")
# The weld group is treated as a collection of lines.
# Origin (0,0) is at the bottom-left corner of the C-shape.
# Weld 1 (top): centroid (25, 100), length 50
# Weld 2 (bottom): centroid (25, 0), length 50
# Weld 3 (vertical): centroid (0, 50), length 100
total_length = w1_length + w2_length + w3_length

# Calculate centroid coordinates (x_G, y_G)
x_G = (w1_length * 25 + w2_length * 25 + w3_length * 0) / total_length
y_G = (w1_length * 100 + w2_length * 0 + w3_length * 50) / total_length

print(f"Total weld length L = {w1_length} + {w2_length} + {w3_length} = {total_length} mm")
print(f"Centroid G = ({x_G:.2f} mm, {y_G:.2f} mm)")
print("-" * 30)

# --- 2. FORCES AND MOMENT ---
print("Step 2: Calculating Forces and Moment")
# The weld throat area A = total_length * t, where t is the throat thickness.
# Primary Shear Stress (acting downwards)
# tau_prime = P / A = P / (total_length * t)
tau_prime_times_t = P / total_length

# The force P is applied at x = 50 + 150 = 200 mm from the origin.
e = (w1_length + eccentricity_arm) - x_G
M = P * e  # Moment about the centroid G

print(f"Direct Shear Force P = {P} N")
print(f"Eccentricity e = ({w1_length} + {eccentricity_arm}) - {x_G:.2f} = {e:.2f} mm")
print(f"Moment M = P * e = {P} * {e:.2f} = {M:.2f} N-mm")
print("-" * 30)


# --- 3. POLAR MOMENT OF INERTIA ---
print("Step 3: Calculating Unit Polar Moment of Inertia (Ju)")
# Using parallel axis theorem: I = I_c + A*d^2
# For a line segment, I_c about axis perpendicular to line is negligible.

# I_xu about the centroidal x-axis (y = y_G)
Ixu_w1 = w1_length * (100 - y_G)**2
Ixu_w2 = w2_length * (0 - y_G)**2
Ixu_w3 = w3_length**3 / 12  # I_c for vertical line about its center
I_xu = Ixu_w1 + Ixu_w2 + Ixu_w3

# I_yu about the centroidal y-axis (x = x_G)
Iyu_w1 = (w1_length**3 / 12) + w1_length * (25 - x_G)**2
Iyu_w2 = (w2_length**3 / 12) + w2_length * (25 - x_G)**2
Iyu_w3 = w3_length * (0 - x_G)**2
I_yu = Iyu_w1 + Iyu_w2 + Iyu_w3

# Unit polar moment of inertia J_u = J / t
J_u = I_xu + I_yu
print(f"Unit Moment of Inertia about x-axis, I_xu = {I_xu:.2f} mm^3")
print(f"Unit Moment of Inertia about y-axis, I_yu = {I_yu:.2f} mm^3")
print(f"Unit Polar Moment of Inertia, J_u = I_xu + I_yu = {J_u:.2f} mm^3")
print("-" * 30)

# --- 4. STRESS ANALYSIS AT CRITICAL POINT ---
print("Step 4: Stress Analysis at Critical Point")
# The critical point is farthest from the centroid. Let's check the top-right corner (50, 100).
crit_x = 50.0
crit_y = 100.0
print(f"Critical point selected at ({crit_x}, {crit_y})")

# r is the vector from the centroid to the critical point
r_x = crit_x - x_G
r_y = crit_y - y_G

# Stresses multiplied by throat thickness 't'
# Primary stress component (downwards)
tau_prime_y_comp = -tau_prime_times_t

# Secondary (torsional) stress components for clockwise moment M
tau_double_prime_x_comp = (M * r_y) / J_u
tau_double_prime_y_comp = -(M * r_x) / J_u

print(f"Primary shear stress component (vertical): τ'_y * t = -{P}/{total_length} = {tau_prime_y_comp:.2f} N/mm")
print(f"Secondary shear stress component (horizontal): τ''_x * t = (M * r_y) / J_u = ({M:.0f} * {r_y:.2f}) / {J_u:.2f} = {tau_double_prime_x_comp:.2f} N/mm")
print(f"Secondary shear stress component (vertical): τ''_y * t = -(M * r_x) / J_u = -({M:.0f} * {r_x:.2f}) / {J_u:.2f} = {tau_double_prime_y_comp:.2f} N/mm")


# Total stress components
tau_total_x_comp = tau_double_prime_x_comp
tau_total_y_comp = tau_prime_y_comp + tau_double_prime_y_comp
print("\nTotal stress components:")
print(f"τ_x * t = {tau_total_x_comp:.2f} N/mm")
print(f"τ_y * t = {tau_prime_y_comp:.2f} + {tau_double_prime_y_comp:.2f} = {tau_total_y_comp:.2f} N/mm")

# Resultant stress
tau_max_times_t = math.sqrt(tau_total_x_comp**2 + tau_total_y_comp**2)

print("\nResultant stress equation:")
print(f"τ_max * t = sqrt((τ_x*t)^2 + (τ_y*t)^2)")
print(f"τ_max * t = sqrt({tau_total_x_comp:.2f}^2 + {tau_total_y_comp:.2f}^2) = {tau_max_times_t:.2f} N/mm")
print("-" * 30)

# --- 5. DETERMINE WELD SIZE ---
print("Step 5: Final Weld Size Calculation")
# Set max stress to permissible stress to find required throat thickness 't'
# tau_max = tau_perm => (tau_max * t) / t = tau_perm
t = tau_max_times_t / tau_perm
print(f"Required throat thickness t = (τ_max * t) / τ_perm = {tau_max_times_t:.2f} / {tau_perm} = {t:.2f} mm")

# For a fillet weld, throat thickness t = 0.707 * s (leg size)
s = t / 0.707
print("\nFinal equation for weld leg size 's':")
print(f"s = t / 0.707 = {t:.2f} / 0.707 = {s:.2f} mm")
print("-" * 30)
print(f"The required weld size (leg) is {s:.2f} mm.")
print("=" * 30)

final_answer = round(s, 2)