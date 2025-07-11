import math

def linear_interpolate(x, x_points, y_points):
    """
    Performs linear interpolation.
    Finds y for a given x, using two points (x1, y1) and (x2, y2).
    x: the point to find the y-value for.
    x_points: a list or tuple of two x-values [x1, x2].
    y_points: a list or tuple of two y-values [y1, y2].
    """
    x1, x2 = x_points
    y1, y2 = y_points
    
    # Handle the case where the slope is vertical (unlikely in this context)
    if x1 == x2:
        return y1 if x == x1 else float('nan')
        
    y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    return y

# --- Given Data ---
# Test conditions
P_plate = 30 * 1000  # Load on plate (N)
d_plate = 305        # Diameter of plate (mm)
w_s_test = 2460 / 1000 # Deflection on subgrade (mm)
w_t_test = 1080 / 1000 # Deflection on trial section (mm)
h_trial = 300        # Thickness of trial pavement (mm)

# Design conditions
W_design = 1.80 * 1000 # Design wheel load mass (kg)
g = 9.81             # Acceleration due to gravity (m/s^2)
p_design_mpa = 600 / 1000 # Tyre pressure (kN/m^2 -> MPa or N/mm^2)
w_allowable = 1.00   # Max allowable deflection (mm)

# --- Step 1: Determine Subgrade Modulus of Elasticity (E_s) ---
print("--- Step 1: Analyzing Subgrade Properties ---")
a_plate = d_plate / 2
p_plate = P_plate / (math.pi * a_plate**2)
# Using the formula for flexible plate deflection: w = (1.5 * p * a) / E
# Rearranging for E: E = (1.5 * p * a) / w
E_s = (1.5 * p_plate * a_plate) / w_s_test
print(f"Plate radius (a_plate) = {a_plate:.2f} mm")
print(f"Plate pressure (p_plate) = {P_plate} N / (pi * {a_plate:.2f}^2 mm^2) = {p_plate:.4f} MPa")
print("Using formula E_s = (1.5 * p_plate * a_plate) / w_s_test:")
print(f"E_s = (1.5 * {p_plate:.4f} MPa * {a_plate:.2f} mm) / {w_s_test:.3f} mm")
print(f"Subgrade Modulus (E_s) = {E_s:.2f} MPa\n")


# --- Step 2: Determine Pavement Modulus (E_p) ---
print("--- Step 2: Analyzing Pavement Properties from Trial Section ---")
# Calculate the deflection factor F2 from the trial section test
F2_test = w_t_test / w_s_test
ha_ratio_test = h_trial / a_plate
print(f"Deflection factor (F2_test) = w_t / w_s = {w_t_test:.3f} mm / {w_s_test:.3f} mm = {F2_test:.4f}")
print(f"Thickness/Radius ratio (h/a_test) = {h_trial} mm / {a_plate:.2f} mm = {ha_ratio_test:.4f}")

# Burmister chart data for F2 vs. K=E_p/E_s at h/a ≈ 2.0
# We will interpolate to find the K that gives F2_test
K_points = [10, 20]
F2_points_at_ha2 = [0.45, 0.32] # F2 values for K=10 and K=20 at h/a=2.0
K_ratio = linear_interpolate(F2_test, F2_points_at_ha2, K_points)
E_p = K_ratio * E_s
print("\nInterpolating from Burmister chart data to find modulus ratio (K = E_p/E_s)...")
print(f"Found Modulus Ratio (K) = {K_ratio:.2f}")
print(f"Pavement Modulus (E_p) = K * E_s = {K_ratio:.2f} * {E_s:.2f} MPa = {E_p:.2f} MPa\n")


# --- Step 3: Analyze Design Load & Required F2 ---
print("--- Step 3: Analyzing Design Load Conditions ---")
P_design = W_design * g
a_design = math.sqrt(P_design / (math.pi * p_design_mpa))
print(f"Design Load (P_design) = {W_design/1000:.2f} ton * {g} m/s^2 = {P_design:.2f} N")
print(f"Design Contact Radius (a_design) = sqrt({P_design:.2f} N / (pi * {p_design_mpa:.3f} N/mm^2)) = {a_design:.2f} mm")

# Calculate hypothetical deflection on subgrade under design load
w_s_design = (1.5 * p_design_mpa * a_design) / E_s
print(f"Hypothetical deflection on subgrade (w_s_design) = (1.5 * {p_design_mpa:.3f} MPa * {a_design:.2f} mm) / {E_s:.2f} MPa = {w_s_design:.3f} mm")

# Calculate the required deflection factor for the design
F2_design = w_allowable / w_s_design
print(f"Required Deflection Factor (F2_design) = w_allowable / w_s_design = {w_allowable:.2f} mm / {w_s_design:.3f} mm = {F2_design:.4f}\n")


# --- Step 4: Determine Required Pavement Thickness (h_design) ---
print("--- Step 4: Determining Required Pavement Thickness ---")
# We need to find the h/a ratio that gives F2_design for our K_ratio.
# This requires 2D interpolation. We'll do it in two steps.
# Burmister data for K=10 and K=20 at two h/a points
ha_points = [1.0, 2.0]
F2_at_K10 = [0.68, 0.45]
F2_at_K20 = [0.48, 0.32]

# First, create an interpolated F2 curve for our specific K_ratio
F2_for_K_ratio_at_ha1 = linear_interpolate(K_ratio, K_points, [F2_at_K10[0], F2_at_K20[0]])
F2_for_K_ratio_at_ha2 = linear_interpolate(K_ratio, K_points, [F2_at_K10[1], F2_at_K20[1]])
print(f"Interpolating to find points on the curve for K = {K_ratio:.2f}:")
print(f"  - At h/a = {ha_points[0]}, F2 = {F2_for_K_ratio_at_ha1:.4f}")
print(f"  - At h/a = {ha_points[1]}, F2 = {F2_for_K_ratio_at_ha2:.4f}")

# Now, interpolate along this new curve to find the h/a for our F2_design
F2_points_for_K = [F2_for_K_ratio_at_ha1, F2_for_K_ratio_at_ha2]
ha_ratio_design = linear_interpolate(F2_design, F2_points_for_K, ha_points)
print(f"\nInterpolating again to find the h/a ratio that gives F2 = {F2_design:.4f}:")
print(f"Required Thickness/Radius Ratio (h/a_design) = {ha_ratio_design:.4f}")

# Final calculation
h_design = ha_ratio_design * a_design
print("\n--- Final Calculation ---")
print("Required Pavement Thickness (h_design) = h/a_design * a_design")
print(f"h_design = {ha_ratio_design:.4f} * {a_design:.2f} mm")
print(f"h_design = {h_design:.2f} mm")

# --- Final Answer ---
final_answer = h_design
print(f"\n<<<The required pavement thickness is {final_answer:.2f} mm.>>>")
print(f"<<<{final_answer:.2f}>>>")