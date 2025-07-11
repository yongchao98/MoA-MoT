import math

# Step 1: Define given data and weld geometry
P = 90 * 1000  # Force in N
tau_perm = 250  # Permissible shear stress in N/mm^2
l1 = 50.0  # Length of weld W1 in mm
l2 = 50.0  # Length of weld W2 in mm
l3 = 100.0 # Length of weld W3 in mm
d = 50.0 + 150.0 # Distance from weld W3 to the force application point

print("--- Step 1: Given Data ---")
print(f"Eccentric Force (P): {P/1000} kN")
print(f"Permissible Shear Stress (τ_perm): {tau_perm} N/mm^2")
print(f"Weld lengths: W1={l1} mm, W2={l2} mm, W3={l3} mm")
print("-" * 30 + "\n")

# Step 2: Calculate the centroid of the weld group
# Origin at the bottom-left corner of the weld group (bottom of W3)
# Centroid of each weld segment:
# c1 for W1, c2 for W2, c3 for W3
x1, y1 = l1 / 2, l3
x2, y2 = l2 / 2, 0
x3, y3 = 0, l3 / 2

L = l1 + l2 + l3 # Total length of the weld

# Centroid coordinates (x_bar, y_bar)
x_bar = (l1 * x1 + l2 * x2 + l3 * x3) / L
y_bar = (l1 * y1 + l2 * y2 + l3 * y3) / L

print("--- Step 2: Centroid of Weld Group (G) ---")
print(f"Total weld length (L): {L} mm")
print(f"Centroid G = (x_bar, y_bar) = ({x_bar:.2f} mm, {y_bar:.2f} mm)")
print("-" * 30 + "\n")

# Step 3: Calculate Moment and Section Properties (per unit throat thickness)
# Eccentricity (e) and Moment (M)
e = d - x_bar
M = P * e

# Calculate unit moment of inertia (J_u = J/t) using parallel axis theorem
# I_u = I_c_u + l*d^2
# For W1 (horizontal):
I_xu1 = (l1 * (y1 - y_bar)**2)
I_yu1 = (l1**3 / 12) + (l1 * (x1 - x_bar)**2)
# For W2 (horizontal):
I_xu2 = (l2 * (y2 - y_bar)**2)
I_yu2 = (l2**3 / 12) + (l2 * (x2 - x_bar)**2)
# For W3 (vertical):
I_xu3 = (l3**3 / 12) + (l3 * (y3 - y_bar)**2)
I_yu3 = (l3 * (x3 - x_bar)**2)

# Total unit moments of inertia
I_xu = I_xu1 + I_xu2 + I_xu3
I_yu = I_yu1 + I_yu2 + I_yu3
J_u = I_xu + I_yu # Unit polar moment of inertia

print("--- Step 3: Moment and Section Properties ---")
print(f"Eccentricity (e): {d} - {x_bar:.2f} = {e:.2f} mm")
print(f"Moment (M = P * e): {M/1e6:.2f} kNm")
print(f"Unit Moment of Inertia about x-axis (I_xu): {I_xu:.2f} mm^3")
print(f"Unit Moment of Inertia about y-axis (I_yu): {I_yu:.2f} mm^3")
print(f"Unit Polar Moment of Inertia (J_u = I_xu + I_yu): {J_u:.2f} mm^3")
print("-" * 30 + "\n")

# Step 4: Calculate Stresses at the most critical point
# The critical points are the corners farthest from the centroid.
# Let's check the bottom-right corner (Point B) of weld W2.
# Coordinates of Point B: (50, 0)
# Coordinates relative to centroid:
r_x = 50.0 - x_bar
r_y = 0.0 - y_bar

# Primary shear stress (τ') acts vertically downwards
# τ' = P / A = P / (L * t). We calculate tau_prime_over_t = P / L
tau_prime_y_over_t = -P / L  # Negative for downwards direction

# Secondary shear stress (τ'') due to torsion
# τ''_x = -M*r_y / J = -(M/J_u) * r_y / t
# τ''_y =  M*r_x / J =  (M/J_u) * r_x / t
tau_double_prime_x_over_t = -M * r_y / J_u
tau_double_prime_y_over_t = M * r_x / J_u

# Total shear stress components (divided by t)
tau_x_over_t = tau_double_prime_x_over_t
tau_y_over_t = tau_prime_y_over_t + tau_double_prime_y_over_t

# Resultant shear stress (divided by t)
K_res = math.sqrt(tau_x_over_t**2 + tau_y_over_t**2)

print("--- Step 4: Stress Analysis at Critical Point B (50, 0) ---")
print(f"Primary Shear Stress (τ') component: τ'_y = {P/1000} kN / ({L} mm * t) = {P/L:.2f}/t N/mm^2")
print(f"Secondary Shear Stress (τ'') components:")
print(f"  τ''_x = -({M:.0f} * {r_y:.2f}) / ({J_u:.0f} * t) = {tau_double_prime_x_over_t:.2f}/t N/mm^2")
print(f"  τ''_y = ({M:.0f} * {r_x:.2f}) / ({J_u:.0f} * t) = {tau_double_prime_y_over_t:.2f}/t N/mm^2")
print(f"Total Resultant Shear Stress (τ_res) = {K_res:.2f} / t N/mm^2")
print("-" * 30 + "\n")

# Step 5: Determine the required weld size (s)
# τ_res = τ_perm => K_res / t = τ_perm
# Solve for throat thickness (t)
t = K_res / tau_perm

# Solve for weld leg size (s), where t = 0.707 * s
s = t / 0.707

print("--- Step 5: Determine Required Weld Size (s) ---")
print("The maximum resultant shear stress must be less than or equal to the permissible stress.")
print("Equation: τ_res = τ_perm")
print(f"{K_res:.2f} / t = {tau_perm}")
print(f"Required throat thickness (t) = {K_res:.2f} / {tau_perm} = {t:.2f} mm")
print("\nFinal equation for weld leg size (s):")
print(f"s = t / 0.707")
print(f"s = {t:.2f} / 0.707")
print(f"s = {s:.2f} mm")
print("-" * 30)

# Final answer
final_answer = round(s, 1)
print(f"\nThe required weld size (leg) is {final_answer} mm.")
print(f'<<<s = {K_res:.1f} / ({tau_perm}*0.707) = {s:.1f}>>>')
<<<14.1>>>