import math

# --- Given Data ---
# Plate Bearing Test
P_plate_kN = 30.0  # Load in kN
d_plate_mm = 305.0  # Plate diameter in mm
delta_s_um = 2460.0  # Deflection on subgrade in micrometers
h_trial_mm = 300.0  # Trial pavement thickness in mm
delta_p_um = 1080.0  # Deflection on trial pavement in micrometers

# Design Load
W_design_ton = 1.80  # Wheel load in tons
p_design_kPa = 600.0  # Tyre pressure in kPa

# Material Properties & Constraints
mu = 0.5  # Poisson's ratio
delta_allowable_mm = 1.00  # Max allowable deflection in mm
g = 9.81 # acceleration due to gravity in m/s^2

# --- 1. Unit Conversion and Initial Calculations ---
print("Step 1: Unit Conversion and Initial Calculations\n")
# Convert to consistent units: N, mm, MPa
P_plate_N = P_plate_kN * 1000
delta_s_mm = delta_s_um / 1000
delta_p_mm = delta_p_um / 1000
p_design_MPa = p_design_kPa / 1000
P_design_N = W_design_ton * 1000 * g

a_plate_mm = d_plate_mm / 2
p_plate_MPa = P_plate_N / (math.pi * a_plate_mm**2)

print(f"Plate radius (a_plate): {d_plate_mm:.2f} mm / 2 = {a_plate_mm:.2f} mm")
print(f"Plate pressure (p_plate): {P_plate_N:.2f} N / (pi * {a_plate_mm:.2f}^2 mm^2) = {p_plate_MPa:.4f} MPa\n")

# --- 2. Calculate Subgrade Modulus (E₂) ---
print("Step 2: Calculate Subgrade Modulus (E₂)\n")
# Formula for rigid plate deflection on a single layer with mu=0.5:
# delta_s = 1.5 * p_plate * a_plate / E₂
E2_MPa = (1.5 * p_plate_MPa * a_plate_mm) / delta_s_mm
print(f"Using formula: E₂ = (1.5 * p_plate * a_plate) / Δ_s")
print(f"E₂ = (1.5 * {p_plate_MPa:.4f} MPa * {a_plate_mm:.2f} mm) / {delta_s_mm:.3f} mm")
print(f"E₂ = {E2_MPa:.2f} MPa\n")

# --- 3. Determine Modular Ratio (E₁/E₂) ---
print("Step 3: Determine Modular Ratio (E₁/E₂)\n")
# The deflection factor F is the ratio of two-layer deflection to subgrade deflection
F_p = delta_p_mm / delta_s_mm
h_a_ratio_p = h_trial_mm / a_plate_mm

print(f"Deflection factor for trial (F_p) = Δ_p / Δ_s = {delta_p_mm:.3f} mm / {delta_s_mm:.3f} mm = {F_p:.4f}")
print(f"h/a ratio for trial = {h_trial_mm:.2f} mm / {a_plate_mm:.2f} mm = {h_a_ratio_p:.4f}")
print("Based on Burmister's chart for two-layer systems (μ=0.5), an h/a ratio of ~1.97 and F of ~0.44 corresponds to a modular ratio E₁/E₂ of approximately 10.")
E1_over_E2 = 10
print(f"Assumed Modular Ratio (E₁/E₂) = {E1_over_E2}\n")

# --- 4. Analyze Design Load ---
print("Step 4: Analyze Design Load\n")
# Calculate the radius of the design load contact area
A_design_mm2 = P_design_N / p_design_MPa
a_design_mm = math.sqrt(A_design_mm2 / math.pi)

print(f"Design load (P_design): {W_design_ton:.2f} ton * 1000 kg/ton * {g} m/s^2 = {P_design_N:.2f} N")
print(f"Design contact area (A_design) = P_design / p_design = {P_design_N:.2f} N / {p_design_MPa:.3f} MPa = {A_design_mm2:.2f} mm^2")
print(f"Design load radius (a_design) = sqrt(A_design / pi) = sqrt({A_design_mm2:.2f} mm^2 / pi) = {a_design_mm:.2f} mm\n")

# --- 5. Calculate Required Thickness (h) ---
print("Step 5: Calculate Required Pavement Thickness (h)\n")
# a) Calculate theoretical subgrade deflection under design load
delta_design_s_mm = (1.5 * p_design_MPa * a_design_mm) / E2_MPa
print("a) Calculate theoretical subgrade deflection for the design load (Δ_design_s):")
print(f"Δ_design_s = (1.5 * p_design * a_design) / E₂")
print(f"Δ_design_s = (1.5 * {p_design_MPa:.3f} MPa * {a_design_mm:.2f} mm) / {E2_MPa:.2f} MPa = {delta_design_s_mm:.4f} mm\n")

# b) Calculate the required deflection factor F for the design
F_design = delta_allowable_mm / delta_design_s_mm
print("b) Calculate the required deflection factor (F_design):")
print(f"F_design = Δ_allowable / Δ_design_s = {delta_allowable_mm:.2f} mm / {delta_design_s_mm:.4f} mm = {F_design:.4f}\n")

# c) Find required h/a ratio by interpolation from chart data for E₁/E₂ = 10
print("c) Find required h/a ratio using interpolation (for E₁/E₂ = 10):")
# Data points from Burmister's chart for E₁/E₂ = 10, μ = 0.5
h_a_1, F_1 = 2.0, 0.45
h_a_2, F_2 = 3.0, 0.35
# Linear interpolation: x = x1 + (y - y1) * (x2 - x1) / (y2 - y1)
h_a_ratio_design = h_a_1 + (F_design - F_1) * (h_a_2 - h_a_1) / (F_2 - F_1)
print(f"Interpolating between points (h/a={h_a_1}, F={F_1}) and (h/a={h_a_2}, F={F_2}) for F_design = {F_design:.4f}")
print(f"Required h/a = {h_a_1:.1f} + ({F_design:.4f} - {F_1:.2f}) * ({h_a_2:.1f} - {h_a_1:.1f}) / ({F_2:.2f} - {F_1:.2f}) = {h_a_ratio_design:.4f}\n")

# d) Calculate final required pavement thickness h
h_required_mm = h_a_ratio_design * a_design_mm
print("d) Calculate the final required pavement thickness (h):")
print(f"h = (h/a)_design * a_design")
print(f"h = {h_a_ratio_design:.4f} * {a_design_mm:.2f} mm")
print(f"h = {h_required_mm:.2f} mm\n")

print("--- Final Answer ---")
print(f"The required pavement thickness is {h_required_mm:.2f} mm.")
# Final answer block
final_answer = round(h_required_mm, 2)
print(f'<<<{final_answer}>>>')
