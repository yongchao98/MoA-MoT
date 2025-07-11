import math

# --- Given Data ---
# Note: All units are converted to Newtons (N) and millimeters (mm) for consistency.

# Plate Bearing Test Parameters
P_test_kN = 30.0  # kN
d_test_mm = 305.0 # mm
delta_subgrade_um = 2460.0 # deflection on subgrade (μm)
delta_pavement_um = 1080.0 # deflection on trial pavement (μm)
h_trial_mm = 300.0 # trial pavement thickness (mm)

# Design Load Parameters
P_design_ton = 1.80 # ton
p_design_kPa = 600.0 # kN/m^2
delta_design_mm = 1.00 # mm

# Material Properties
poisson_ratio = 0.5 # μ

# --- Helper Functions for Interpolation (Simulating Burmister Chart Lookup) ---

def interpolate(x, x1, y1, x2, y2):
    """Linearly interpolates to find y for a given x."""
    return y1 + (x - x1) * (y2 - y1) / (x2 - x1)

def inv_interpolate(y, x1, y1, x2, y2):
    """Linearly interpolates to find x for a given y."""
    # To avoid division by zero if y1 == y2
    if (y2 - y1) == 0:
        return x1
    return x1 + (y - y1) * (x2 - x1) / (y2 - y1)

# --- Calculations ---

print("### Step 1: Determine Subgrade Modulus (E₂) ###")

# Convert units
P_test_N = P_test_kN * 1000
delta_subgrade_mm = delta_subgrade_um / 1000

# Calculate test plate radius and pressure
a_test_mm = d_test_mm / 2
area_test = math.pi * a_test_mm**2
p_test_MPa = P_test_N / area_test

# Calculate E₂ using the single-layer deflection formula: delta = 1.5 * p * a / E₂
E2_MPa = (1.5 * p_test_MPa * a_test_mm) / delta_subgrade_mm

print(f"Test plate radius (a_test): {a_test_mm:.2f} mm")
print(f"Test pressure (p_test): {p_test_MPa:.4f} N/mm² (MPa)")
print("E₂ = (1.5 * p_test * a_test) / Δ_subgrade")
print(f"E₂ = (1.5 * {p_test_MPa:.4f} * {a_test_mm:.2f}) / {delta_subgrade_mm:.4f}")
print(f"Subgrade Modulus (E₂): {E2_MPa:.2f} N/mm² (MPa)\n")


print("### Step 2: Determine Pavement Modulus (E₁) ###")

# Convert trial section deflection to mm
delta_pavement_mm = delta_pavement_um / 1000

# Calculate the deflection factor F₂ for the trial section
# Formula: delta = (1.5 * p * a / E₂) * F₂  =>  F₂ = (delta * E₂) / (1.5 * p * a)
F2_trial = (delta_pavement_mm * E2_MPa) / (1.5 * p_test_MPa * a_test_mm)
h_a_ratio_trial = h_trial_mm / a_test_mm

print(f"Trial h/a ratio (h_trial / a_test): {h_a_ratio_trial:.3f}")
print("F₂_trial = (Δ_pavement * E₂) / (1.5 * p_test * a_test)")
print(f"F₂_trial = ({delta_pavement_mm:.4f} * {E2_MPa:.2f}) / (1.5 * {p_test_MPa:.4f} * {a_test_mm:.2f})")
print(f"Calculated Deflection Factor (F₂_trial): {F2_trial:.4f}\n")

# Interpolate from standard Burmister charts to find E₁/E₂
# We use data points from the chart for h/a ≈ 2.0 (since our h/a is 1.967)
# Chart data at h/a=2: Point1(E₁/E₂=10, F₂=0.48), Point2(E₁/E₂=20, F₂=0.38)
E1_E2_ratio = inv_interpolate(F2_trial, 10, 0.48, 20, 0.38)
E1_MPa = E1_E2_ratio * E2_MPa

print("Using interpolation from standard Burmister chart data to find E₁/E₂ ratio...")
print(f"Found E₁/E₂ Ratio: {E1_E2_ratio:.2f}")
print(f"Pavement Modulus (E₁): {E1_E2_ratio:.2f} * {E2_MPa:.2f} = {E1_MPa:.2f} N/mm² (MPa)\n")


print("### Step 3: Determine Design Load Parameters ###")

# Convert design load and pressure to N and mm
P_design_N = P_design_ton * 1000 * 9.81  # using g = 9.81 m/s²
p_design_MPa = p_design_kPa / 1000 # 600 kN/m^2 = 0.6 N/mm^2 (MPa)

# Calculate design contact radius from Area = Load / Pressure
a_design_mm = math.sqrt(P_design_N / (p_design_MPa * math.pi))

print(f"Design wheel load: {P_design_N:.2f} N")
print(f"Design tyre pressure: {p_design_MPa:.2f} N/mm² (MPa)")
print("a_design = sqrt(P_design / (p_design * π))")
print(f"a_design = sqrt({P_design_N:.2f} / ({p_design_MPa:.2f} * π))")
print(f"Design Contact Radius (a_design): {a_design_mm:.2f} mm\n")


print("### Step 4 & 5: Find Required Thickness-to-Radius Ratio (h/a) ###")
# Calculate the required deflection factor F₂ for the design scenario
F2_design = (delta_design_mm * E2_MPa) / (1.5 * p_design_MPa * a_design_mm)
print("F₂_design = (Δ_design * E₂) / (1.5 * p_design * a_design)")
print(f"F₂_design = ({delta_design_mm:.2f} * {E2_MPa:.2f}) / (1.5 * {p_design_MPa:.2f} * {a_design_mm:.2f})")
print(f"Required Deflection Factor (F₂_design): {F2_design:.4f}\n")

# Now find the h/a ratio that gives this F₂ for our E₁/E₂ ratio.
# We need to interpolate along the E₁/E₂ = 14.1 (our calculated value) curve.
# First, find points on that curve by interpolating known chart data.
# Chart data for E₁/E₂=10: (h/a=1, F₂=0.72), (h/a=2, F₂=0.48)
# Chart data for E₁/E₂=20: (h/a=1, F₂=0.60), (h/a=2, F₂=0.38)
F2_at_ha1 = interpolate(E1_E2_ratio, 10, 0.72, 20, 0.60)
F2_at_ha2 = interpolate(E1_E2_ratio, 10, 0.48, 20, 0.38)
print(f"Interpolating to find F₂ points for the curve E₁/E₂ = {E1_E2_ratio:.2f}:")
print(f"- At h/a = 1.0, F₂ ≈ {F2_at_ha1:.4f}")
print(f"- At h/a = 2.0, F₂ ≈ {F2_at_ha2:.4f}\n")


# Now find the h/a that corresponds to F₂_design using our two new points
h_a_ratio_design = inv_interpolate(F2_design, 1, F2_at_ha1, 2, F2_at_ha2)
print("Using these points, we find the design h/a ratio for our required F₂:")
print(f"Required Design h/a Ratio: {h_a_ratio_design:.3f}\n")


print("### Step 6: Calculate Final Required Pavement Thickness (h_design) ###")
h_design_mm = h_a_ratio_design * a_design_mm
print("h_design = (h/a)_design * a_design")
print(f"h_design = {h_a_ratio_design:.3f} * {a_design_mm:.2f}")
print(f"Final Required Pavement Thickness: {h_design_mm:.1f} mm")

# --- Final Answer ---
final_answer = round(h_design_mm, 1)
print(f"\n<<< {final_answer} >>>")