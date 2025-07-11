import math

# --- Helper function for linear interpolation ---
def interpolate(x, x1, y1, x2, y2):
    """Linearly interpolates to find y for a given x."""
    if x2 == x1:
        return y1
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1)

# --- Given Data ---
# Test 1 (on subgrade)
P1 = 30e3  # Load (N)
d1 = 305e-3 # Plate diameter (m)
a1 = d1 / 2 # Plate radius (m)
delta1 = 2460e-6 # Deflection (m)

# Test 2 (on pavement)
P2 = 30e3 # Load (N)
h_trial = 300e-3 # Trial pavement thickness (m)
delta2 = 1080e-6 # Deflection (m)

# Design Parameters
W_design = 1.80 * 1000 # Design wheel load weight (kg)
P_design = W_design * 9.81 # Design wheel load (N)
q_design = 600e3 # Tyre pressure (N/m^2 or Pa)
delta_design_limit = 1.00e-3 # Max allowable deflection (m)

# Material Property
mu = 0.5 # Poisson's ratio

# --- Burmister F2 Chart Data (for flexible load, mu=0.5) ---
# F2_chart[E1/E2][h/a] = F2
F2_chart = {
    5: {1: 0.78, 2: 0.53, 3: 0.40, 4: 0.31},
    10: {1: 0.61, 2: 0.42, 3: 0.30, 4: 0.23},
    20: {1: 0.47, 2: 0.31, 3: 0.22, 4: 0.17}
}

# --- Step 1: Calculate Subgrade Modulus (E2) ---
# Using rigid plate formula: delta = 1.18 * q * a / E
# q = P / (pi * a^2)
q1 = P1 / (math.pi * a1**2)
E2 = (1.18 * q1 * a1) / delta1
print(f"Step 1: Calculated Subgrade Modulus (E2) = {E2/1e6:.2f} MPa")
print("-" * 20)

# --- Step 2: Determine E1/E2 Ratio ---
# Deflection factor F2 from test data
# Note: F2 is the ratio of two-layer deflection to single-layer deflection
F2_trial = delta2 / delta1
h_a_trial = h_trial / a1
print(f"Step 2: From trial data...")
print(f"Deflection Factor (F2) = {F2_trial:.4f}")
print(f"h/a ratio = {h_a_trial:.4f}")

# Interpolate to find E1/E2 ratio from the F2 chart.
# We interpolate along the h/a=2 line, as h_a_trial=1.967 is very close to 2.
h_a_for_interp = 2
e_ratios = sorted(F2_chart.keys()) # [5, 10, 20]
f2_values_at_h_a = [F2_chart[ratio][h_a_for_interp] for ratio in e_ratios]

# Find which two points F2_trial lies between
e1_ratio, e2_ratio = 0, 0
f2_1, f2_2 = 0, 0
for i in range(len(e_ratios) - 1):
    # F2 values decrease as E1/E2 increases
    if f2_values_at_h_a[i+1] <= F2_trial <= f2_values_at_h_a[i]:
        e1_ratio, e2_ratio = e_ratios[i], e_ratios[i+1]
        f2_1, f2_2 = f2_values_at_h_a[i], f2_values_at_h_a[i+1]
        break

# Interpolate: y is E1/E2, x is F2
E1_E2_ratio = interpolate(F2_trial, f2_1, e1_ratio, f2_2, e2_ratio)
print(f"Calculated E1/E2 Ratio = {E1_E2_ratio:.4f}")
print("-" * 20)

# --- Step 3: Calculate Design Load Radius (a) ---
area_design = P_design / q_design
a_design = math.sqrt(area_design / math.pi)
print(f"Step 3: Design load contact radius (a) = {a_design * 1000:.2f} mm")
print("-" * 20)

# --- Step 4: Determine Required Pavement Thickness (h) ---
# For flexible tyre load, delta = 1.5 * q * a / E2 * F2
# We need to find the F2 required to meet the design deflection limit.
F2_design_req = (delta_design_limit * E2) / (1.5 * q_design * a_design)
print(f"Step 4: Solving for required thickness (h)...")
print(f"Required Deflection Factor (F2) for design = {F2_design_req:.4f}")

# Now, find h/a for our E1_E2_ratio that gives F2_design_req.
# We'll create a new set of (h/a, F2) points for our specific E1_E2_ratio.
h_a_points = sorted(F2_chart[e1_ratio].keys()) # [1, 2, 3, 4]
F2_points_for_E1_E2 = []
for ha_point in h_a_points:
    f2_at_e1 = F2_chart[e1_ratio][ha_point]
    f2_at_e2 = F2_chart[e2_ratio][ha_point]
    # Interpolate: y is F2, x is E1/E2
    f2_interp = interpolate(E1_E2_ratio, e1_ratio, f2_at_e1, e2_ratio, f2_at_e2)
    F2_points_for_E1_E2.append(f2_interp)

# Find which two points F2_design_req lies between
ha1, ha2 = 0, 0
f2_at_ha1, f2_at_ha2 = 0, 0
for i in range(len(h_a_points) - 1):
    # F2 values decrease as h/a increases
    if F2_points_for_E1_E2[i+1] <= F2_design_req <= F2_points_for_E1_E2[i]:
        ha1, ha2 = h_a_points[i], h_a_points[i+1]
        f2_at_ha1 = F2_points_for_E1_E2[i]
        f2_at_ha2 = F2_points_for_E1_E2[i+1]
        break

# Interpolate: y is h/a, x is F2
h_a_final = interpolate(F2_design_req, f2_at_ha1, ha1, f2_at_ha2, ha2)
print(f"Resulting h/a ratio = {h_a_final:.4f}")

# Final calculation for h
h_final_m = h_a_final * a_design
h_final_mm = h_final_m * 1000

print("-" * 20)
print("Final Calculation:")
print("Required pavement thickness (h) = (h/a ratio) * (design load radius)")
print(f"h = {h_a_final:.4f} * {a_design*1000:.2f} mm")
print(f"h = {h_final_mm:.1f} mm")
print("-" * 20)

# Final answer in specified format
final_answer_mm = round(h_final_mm, 1)
print(f'<<<An equation showing the final calculation is h = {h_a_final:.3f} * {a_design*1000:.2f} mm = {h_final_mm:.1f} mm>>>')
print(f'<<<{final_answer_mm}>>>')
