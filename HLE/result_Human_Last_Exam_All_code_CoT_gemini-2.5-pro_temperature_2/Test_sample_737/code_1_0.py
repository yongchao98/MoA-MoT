import math

# --- Given Parameters ---
# Plate Bearing Test on Subgrade
delta_subgrade_um = 2460.0
P_plate_kN = 30.0
d_plate_mm = 305.0

# Plate Bearing Test on Trial Section
delta_trial_um = 1080.0
h_trial_mm = 300.0

# Design Parameters
W_design_ton = 1.80
p_design_kPa = 600.0
delta_allowable_mm = 1.00
mu = 0.5
g = 9.81 # m/s^2, standard gravity

# --- Unit Conversions ---
delta_subgrade_mm = delta_subgrade_um / 1000.0
P_plate_N = P_plate_kN * 1000.0
a_plate_mm = d_plate_mm / 2.0
delta_trial_mm = delta_trial_um / 1000.0
P_design_N = W_design_ton * 1000.0 * g
p_design_MPa = p_design_kPa / 1000.0 # 1 MPa = 1 N/mm^2

print("Step-by-step calculation to determine the required pavement thickness:")
print("-" * 70)

# --- Step 1: Determine the subgrade modulus of elasticity (E2) ---
# Formula for surface deflection on a single-layer system (flexible plate, mu=0.5):
# Δ = 1.5 * P / (π * a * E)
# Rearranged for E: E2 = 1.5 * P_plate / (π * a_plate * Δ_subgrade)
E2 = (1.5 * P_plate_N) / (math.pi * a_plate_mm * delta_subgrade_mm)
print("1. Calculate the subgrade modulus of elasticity (E2):")
print(f"E2 = (1.5 * P_plate) / (π * a_plate * Δ_subgrade)")
print(f"E2 = (1.5 * {P_plate_N:.2f} N) / (π * {a_plate_mm:.2f} mm * {delta_subgrade_mm:.3f} mm)")
print(f"E2 = {E2:.2f} N/mm² (MPa)\n")

# --- Step 2: Determine the parameters from the trial section test ---
# The two-layer deflection factor F relates the deflection of the two-layer system (Δ_trial)
# to the deflection of the subgrade alone (Δ_subgrade) for the same load and plate.
# F_trial = Δ_trial / Δ_subgrade
F_trial = delta_trial_mm / delta_subgrade_mm
print("2. Calculate the deflection factor (F_trial) from the trial section test:")
print(f"F_trial = Δ_trial / Δ_subgrade")
print(f"F_trial = {delta_trial_mm:.3f} mm / {delta_subgrade_mm:.3f} mm")
print(f"F_trial = {F_trial:.4f}\n")

# Calculate the h/a ratio for the trial test
h_a_ratio_trial = h_trial_mm / a_plate_mm
print("3. Calculate the h/a ratio for the trial section:")
print(f"h/a (trial) = h_trial / a_plate")
print(f"h/a (trial) = {h_trial_mm:.2f} mm / {a_plate_mm:.2f} mm")
print(f"h/a (trial) = {h_a_ratio_trial:.4f}\n")


# --- Step 3: Determine the parameters for the design load ---
# Calculate the contact radius 'a' for the design wheel load.
# Pressure p = P / (π * a^2), so a = sqrt(P / (p * π))
a_design_mm = math.sqrt(P_design_N / (p_design_MPa * math.pi))
print("4. Calculate the contact radius (a_design) for the design wheel load:")
print(f"a_design = sqrt((Weight * g) / (p_design * π))")
print(f"a_design = sqrt(({W_design_ton * 1000:.2f} kg * {g}) / ({p_design_MPa:.3f} N/mm² * π))")
print(f"a_design = sqrt({P_design_N:.2f} N / {p_design_MPa * math.pi:.4f} N/mm²)")
print(f"a_design = {a_design_mm:.2f} mm\n")

# --- Step 4: Determine the required deflection factor for the design ---
# Using the two-layer deflection formula for the design case:
# Δ_allowable = 1.5 * p_design * a_design / E2 * F_design
# Rearranged for F_design: F_design = (Δ_allowable * E2) / (1.5 * p_design * a_design)
F_design = (delta_allowable_mm * E2) / (1.5 * p_design_MPa * a_design_mm)
print("5. Calculate the required deflection factor (F_design):")
print(f"F_design = (Δ_allowable * E2) / (1.5 * p_design * a_design)")
print(f"F_design = ({delta_allowable_mm:.2f} mm * {E2:.2f} N/mm²) / (1.5 * {p_design_MPa:.3f} N/mm² * {a_design_mm:.2f} mm)")
print(f"F_design = {F_design:.4f}\n")


# --- Step 5: Determine the required pavement thickness (h_design) ---
# The deflection factor F is a function of E1/E2 and h/a.
# Since F_trial ({F_trial:.4f}) is approximately equal to F_design ({F_design:.4f}), we can assume their
# corresponding h/a ratios are also approximately equal.
# h_design / a_design ≈ h_trial / a_plate
# Rearranged for h_design: h_design = a_design * (h_trial / a_plate)
h_design_mm = a_design_mm * h_a_ratio_trial
print("6. Calculate the required pavement thickness (h_design):")
print("Since F_trial ≈ F_design, we can assume (h_design / a_design) ≈ (h_trial / a_plate).")
print(f"h_design = a_design * (h_trial / a_plate)")
print(f"h_design = {a_design_mm:.2f} mm * ({h_trial_mm:.2f} mm / {a_plate_mm:.2f} mm)")
print(f"h_design = {a_design_mm:.2f} mm * {h_a_ratio_trial:.4f}")

print("-" * 70)
print(f"The final required pavement thickness is {h_design_mm:.2f} mm.")
print(f"<<<{h_design_mm:.2f}>>>")