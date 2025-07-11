import math

# --- 1. Initial Parameters ---
P = 90 * 1000  # Force in N
tau_perm = 250  # Permissible shear stress in N/mm^2
l1 = 50.0  # Length of top weld W1 in mm
l2 = 50.0  # Length of bottom weld W2 in mm
l3 = 100.0 # Length of vertical weld W3 in mm
e_dist = 150.0 # Distance from weld end to force in mm

print("--- Step 1: Define Problem Parameters ---")
print(f"Force P = {P/1000} kN = {P} N")
print(f"Permissible shear stress τ_perm = {tau_perm} N/mm^2")
print(f"Weld lengths: l1 = {l1} mm, l2 = {l2} mm, l3 = {l3} mm")
print("-" * 40)

# --- 2. Geometric Properties of the Weld Group ---
# Let the origin (0,0) be at the bottom of the vertical weld W3.
# Centroids of individual welds:
# W1: (l1/2, l3) = (25, 100)
# W2: (l2/2, 0) = (25, 0)
# W3: (0, l3/2) = (0, 50)

# Total length of weld
L = l1 + l2 + l3

# Centroid of the weld group (x_bar, y_bar)
x_bar = (l1 * (l1 / 2) + l2 * (l2 / 2) + l3 * 0) / L
y_bar = (l1 * l3 + l2 * 0 + l3 * (l3 / 2)) / L

print("--- Step 2: Calculate Weld Group Centroid (G) ---")
print(f"Total weld length L = {l1} + {l2} + {l3} = {L} mm")
print(f"x_bar = ({l1}*{l1/2} + {l2}*{l2/2} + {l3}*0) / {L} = {x_bar:.2f} mm")
print(f"y_bar = ({l1}*{l3} + {l2}*0 + {l3}*{l3/2}) / {L} = {y_bar:.2f} mm")
print(f"Centroid G is at ({x_bar:.2f}, {y_bar:.2f}) mm")
print("-" * 40)

# --- 3. Calculate Moment of Inertia (J) ---
# We use the parallel axis theorem: I = I_c + A*d^2
# We treat welds as lines, so I_c is (1/12)*t*l^3 for axis perp. to length
# and negligible for axis parallel to length. A = l*t.
# I_xx and I_yy are per unit throat thickness 't'.

# I_xx about the centroid G
I_xx1 = (l1 * (l3 - y_bar)**2)  # For W1
I_xx2 = (l2 * (0 - y_bar)**2)   # For W2
I_xx3 = (1/12 * l3**3) + l3 * (l3/2 - y_bar)**2 # For W3
I_xx = I_xx1 + I_xx2 + I_xx3

# I_yy about the centroid G
I_yy1 = (1/12 * l1**3) + l1 * (l1/2 - x_bar)**2 # For W1
I_yy2 = (1/12 * l2**3) + l2 * (l2/2 - x_bar)**2 # For W2
I_yy3 = (l3 * (0 - x_bar)**2) # For W3
I_yy = I_yy1 + I_yy2 + I_yy3

# Polar moment of inertia J (per unit throat thickness 't')
J_per_t = I_xx + I_yy

print("--- Step 3: Calculate Polar Moment of Inertia (J) ---")
print(f"I_xx (per t) = {I_xx:.2f} mm^3")
print(f"I_yy (per t) = {I_yy:.2f} mm^3")
print(f"J (per t) = I_xx + I_yy = {I_xx:.2f} + {I_yy:.2f} = {J_per_t:.2f} mm^3")
print("-" * 40)

# --- 4. Calculate Stresses ---
# Eccentricity e
e = (l1 + e_dist) - x_bar
# Moment M
M = P * e

# Primary shear stress (per unit throat thickness 't')
tau_prime_per_t = P / L # This is τ' * t

print("--- Step 4: Calculate Forces and Stresses ---")
print("The force P causes a direct shear and a moment about the centroid G.")
print(f"Eccentricity e = ({l1} + {e_dist}) - {x_bar:.2f} = {e:.2f} mm")
print(f"Moment M = P * e = {P} N * {e:.2f} mm = {M:.2f} N-mm")
print(f"Primary shear stress τ' = P / A = {P} / ({L}*t) = {tau_prime_per_t:.2f} / t N/mm^2")
print("\nThe maximum stress occurs at the point farthest from the centroid.")
print("This is at the top-right or bottom-right corner of the weld.")
# Let's analyze the top-right corner: (l1, l3) = (50, 100)
# Coordinates relative to centroid G:
r_x = l1 - x_bar
r_y = l3 - y_bar
print(f"Critical point coordinates relative to G: r_x = {r_x:.2f} mm, r_y = {r_y:.2f} mm")

# Secondary shear stress components (per unit throat thickness 't')
tau_sec_x_per_t = -M * r_y / J_per_t
tau_sec_y_per_t = M * r_x / J_per_t

print("\nSecondary shear stress (τ'') components:")
print(f"τ''_x = -M*r_y / J = -({M:.2f} * {r_y:.2f}) / ({J_per_t:.2f}*t) = {tau_sec_x_per_t:.2f} / t N/mm^2")
print(f"τ''_y =  M*r_x / J =  ({M:.2f} * {r_x:.2f}) / ({J_per_t:.2f}*t) = {tau_sec_y_per_t:.2f} / t N/mm^2")

# Total shear stress components (per unit throat thickness 't')
# Primary shear acts downwards (-y direction)
total_tau_x_per_t = tau_sec_x_per_t
total_tau_y_per_t = tau_sec_y_per_t - tau_prime_per_t

print("\nTotal shear stress components:")
print(f"τ_x = τ''_x = {total_tau_x_per_t:.2f} / t N/mm^2")
print(f"τ_y = τ''_y + τ' = ({tau_sec_y_per_t:.2f}/t) - ({tau_prime_per_t:.2f}/t) = {total_tau_y_per_t:.2f} / t N/mm^2")
print("-" * 40)

# --- 5. Determine Weld Size ---
# Resultant shear stress (per unit throat thickness 't')
tau_max_per_t = math.sqrt(total_tau_x_per_t**2 + total_tau_y_per_t**2)

# Required throat thickness t
t = tau_max_per_t / tau_perm

# Required leg size s
s = t / 0.707

print("--- Step 5: Determine Required Weld Size ---")
print("The maximum resultant shear stress is the vector sum of the components.")
print(f"τ_max = sqrt( (τ_x)^2 + (τ_y)^2 ) = sqrt( ({total_tau_x_per_t:.2f}/t)^2 + ({total_tau_y_per_t:.2f}/t)^2 )")
print(f"τ_max = {tau_max_per_t:.2f} / t N/mm^2")
print("\nSet τ_max = τ_perm to find the required throat thickness 't'.")
print(f"t = {tau_max_per_t:.2f} / τ_perm = {tau_max_per_t:.2f} / {tau_perm} = {t:.3f} mm")
print("\nFinally, calculate the weld leg size 's'.")
print(f"s = t / 0.707 = {t:.3f} / 0.707 = {s:.3f} mm")
print("-" * 40)
print(f"The required size of the weld (leg) is {s:.1f} mm.")

# Final answer in the requested format
print(f"\n<<<{s:.1f}>>>")