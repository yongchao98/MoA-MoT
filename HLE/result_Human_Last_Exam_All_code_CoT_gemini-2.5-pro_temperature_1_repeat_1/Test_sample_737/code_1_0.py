import math

def linear_interp(x_target, x_points, y_points):
    """
    Performs linear interpolation for a target x value given a set of (x, y) points.
    Assumes that the x_points list is monotonic (either increasing or decreasing).
    """
    # Handle cases where x_points are in descending order
    if x_points[0] > x_points[1]:
        # Reverse both lists to make x_points ascending
        x_points = x_points[::-1]
        y_points = y_points[::-1]

    # Handle edge cases
    if x_target <= x_points[0]:
        return y_points[0]
    if x_target >= x_points[-1]:
        return y_points[-1]

    # Find the two points that bracket the target
    for i in range(len(x_points) - 1):
        if x_points[i] <= x_target < x_points[i+1]:
            x1, y1 = x_points[i], y_points[i]
            x2, y2 = x_points[i+1], y_points[i+1]
            
            # Perform linear interpolation
            return y1 + (x_target - x1) * (y2 - y1) / (x2 - x1)
    
    # Fallback in case of unexpected input
    return y_points[-1]

# --- Input Parameters and Constants ---
# Gravitational acceleration
g = 9.81  # m/s^2

# Test 1: Plate bearing test on subgrade
P_test_N = 30 * 1000  # Load in Newtons
d_plate_mm = 305      # Plate diameter in mm
a_plate_mm = d_plate_mm / 2 # Plate radius in mm
delta_subgrade_mm = 2460 / 1000 # Deflection in mm

# Test 2: Plate bearing test on trial pavement
h_trial_mm = 300      # Trial pavement thickness in mm
delta_trial_mm = 1080 / 1000 # Deflection in mm

# Design Parameters
W_design_ton = 1.80   # Design wheel load in tons
p_design_kpa = 600    # Design tyre pressure in kPa
delta_allowable_mm = 1.00 # Allowable deflection in mm

print("### Pavement Thickness Calculation using Burmister's Two-Layer Theory ###\n")

# --- Step 1: Determine Subgrade Modulus (Es) ---
print("--- Step 1: Determine Subgrade Modulus (Es) from Plate Bearing Test ---")
p_plate_mpa = P_test_N / (math.pi * a_plate_mm**2)
print(f"Plate radius, a = {a_plate_mm:.2f} mm")
print(f"Applied pressure under plate, p = {p_plate_mpa:.4f} MPa")
print(f"Measured deflection on subgrade, Δ_subgrade = {delta_subgrade_mm:.2f} mm")

# From Δ = 1.5 * p * a / Es (since F2 = 1 for a single layer)
Es_mpa = (1.5 * p_plate_mpa * a_plate_mm) / delta_subgrade_mm
print(f"Calculated Subgrade Modulus, Es = (1.5 * {p_plate_mpa:.4f} MPa * {a_plate_mm:.2f} mm) / {delta_subgrade_mm:.2f} mm = {Es_mpa:.2f} MPa\n")

# --- Step 2: Determine Pavement to Subgrade Modulus Ratio (Ep/Es) ---
print("--- Step 2: Determine Modulus Ratio (Ep/Es) from Trial Section Test ---")
print(f"Measured deflection on trial section, Δ_trial = {delta_trial_mm:.2f} mm")
# F2 = Δ_trial / Δ_subgrade
F2_trial = delta_trial_mm / delta_subgrade_mm
print(f"Deflection Factor, F2_trial = Δ_trial / Δ_subgrade = {delta_trial_mm:.2f} / {delta_subgrade_mm:.2f} = {F2_trial:.4f}")

h_a_ratio_trial = h_trial_mm / a_plate_mm
print(f"Thickness-to-radius ratio for trial, h/a = {h_trial_mm:.1f} / {a_plate_mm:.2f} = {h_a_ratio_trial:.4f}")

# Interpolate from Burmister chart data to find Ep/Es
print("Interpolating from Burmister's F2 chart data...")
h_a_chart_pts = [1.5, 2.0]
Ep_Es_chart_pts = [10, 20]
F2_chart_data = {
    10: [0.51, 0.47],  # F2 for Ep/Es=10 at h/a=1.5 and 2.0
    20: [0.43, 0.39]   # F2 for Ep/Es=20 at h/a=1.5 and 2.0
}
f2_at_ha_for_epes10 = linear_interp(h_a_ratio_trial, h_a_chart_pts, F2_chart_data[10])
f2_at_ha_for_epes20 = linear_interp(h_a_ratio_trial, h_a_chart_pts, F2_chart_data[20])

# Interpolate Ep/Es based on the target F2_trial (using log scale for Ep/Es)
log_Ep_Es_chart_pts = [math.log10(p) for p in Ep_Es_chart_pts]
f2_pts_for_interp = [f2_at_ha_for_epes20, f2_at_ha_for_epes10] # Reversed order as F2 decreases with Ep/Es
log_Ep_Es_ratio = linear_interp(F2_trial, f2_pts_for_interp, log_Ep_Es_chart_pts)
Ep_Es_ratio = 10**log_Ep_Es_ratio
print(f"For a target F2 of {F2_trial:.4f}, the interpolated modulus ratio, Ep/Es = {Ep_Es_ratio:.2f}\n")

# --- Step 3: Determine Required Pavement Thickness (h_design) ---
print("--- Step 3: Determine Required Pavement Thickness (h_design) ---")
P_design_N = W_design_ton * 1000 * g
p_design_mpa = p_design_kpa / 1000
a_design_mm = math.sqrt(P_design_N / (math.pi * p_design_mpa))
print(f"Design wheel load, P_design = {P_design_N/1000:.2f} kN")
print(f"Design tyre pressure, p_design = {p_design_mpa:.2f} MPa")
print(f"Radius of contact area, a_design = sqrt({P_design_N:.2f} N / (π * {p_design_mpa:.2f} MPa)) = {a_design_mm:.2f} mm")

# Find the required F2 for the design case
F2_design = (delta_allowable_mm * Es_mpa) / (1.5 * p_design_mpa * a_design_mm)
print(f"Allowable deflection, Δ_allowable = {delta_allowable_mm:.2f} mm")
print(f"Required Deflection Factor, F2_design = ({delta_allowable_mm:.2f} * {Es_mpa:.2f}) / (1.5 * {p_design_mpa:.2f} * {a_design_mm:.2f}) = {F2_design:.4f}")

# Interpolate from chart data again to find h/a ratio for the design
h_a_full_chart_pts = [1.0, 1.5, 2.0, 3.0, 4.0]
F2_full_chart_data = {
    10: [0.58, 0.51, 0.47, 0.41, 0.37],
    20: [0.50, 0.43, 0.39, 0.33, 0.30]
}

# For the calculated Ep_Es_ratio, create an F2 curve vs h/a by interpolating
f2_curve_for_design_epes = []
for i in range(len(h_a_full_chart_pts)):
    f2_vals_at_ha = [F2_full_chart_data[20][i], F2_full_chart_data[10][i]]
    f2_interp = linear_interp(math.log10(Ep_Es_ratio), log_Ep_Es_chart_pts, f2_vals_at_ha)
    f2_curve_for_design_epes.append(f2_interp)

# Interpolate h/a from our new curve using F2_design
h_a_design_ratio = linear_interp(F2_design, f2_curve_for_design_epes, h_a_full_chart_pts)
print(f"For Ep/Es = {Ep_Es_ratio:.2f}, the required h/a ratio for F2 = {F2_design:.4f} is {h_a_design_ratio:.4f}")

# --- Final Result ---
h_design_mm = h_a_design_ratio * a_design_mm
print("\n--- Final Required Pavement Thickness ---")
print(f"h_design = h/a_design * a_design")
print(f"h_design = {h_a_design_ratio:.4f} * {a_design_mm:.2f} mm")
print(f"h_design = {h_design_mm:.2f} mm")

print(f"\n<<<{h_design_mm:.2f}>>>")