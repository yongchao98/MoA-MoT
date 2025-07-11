import math

# --- Given Data ---

# Plate Bearing Test (on Subgrade)
load_plate_kN = 30.0
diameter_plate_mm = 305.0
deflection_subgrade_um = 2460.0

# Plate Bearing Test (on Trial Pavement)
h_trial_mm = 300.0
deflection_trial_um = 1080.0

# Design Parameters
design_load_ton = 1.80
tyre_pressure_kPa = 600.0
max_deflection_mm = 1.00
g = 9.81 # m/s^2

# --- Calculations ---

# Convert units to be consistent (N, mm, MPa)
load_plate_N = load_plate_kN * 1000
deflection_subgrade_mm = deflection_subgrade_um / 1000
a_plate_mm = diameter_plate_mm / 2
deflection_trial_mm = deflection_trial_um / 1000

# Step 1: Calculate the Subgrade Modulus (E2)
# Using the formula for deflection on a single layer (elastic half-space):
# Δ = 1.5 * P / (π * a * E)  (for μ=0.5)
# Rearranging for E2: E2 = 1.5 * P / (π * a * Δ)
E2_MPa = (1.5 * load_plate_N) / (math.pi * a_plate_mm * deflection_subgrade_mm)

print("--- Step 1: Calculate Subgrade Modulus (E2) ---")
print(f"The radius of the loading plate is {a_plate_mm:.1f} mm.")
print("Using the formula E2 = (1.5 * P) / (π * a * Δ_subgrade):")
print(f"E2 = (1.5 * {load_plate_N:.0f} N) / (π * {a_plate_mm:.1f} mm * {deflection_subgrade_mm:.3f} mm)")
print(f"E2 = {E2_MPa:.2f} MPa\n")

# Step 2: Calculate Design Load Parameters
load_design_N = design_load_ton * 1000 * g  # Convert tons (mass) to Newtons (force)
pressure_design_MPa = tyre_pressure_kPa / 1000 # Convert kPa to MPa (N/mm^2)

# Area = Force / Pressure
area_design_mm2 = load_design_N / pressure_design_MPa
# Area = π * a^2, so a = sqrt(Area / π)
a_design_mm = math.sqrt(area_design_mm2 / math.pi)

print("--- Step 2: Calculate Design Load Radius (a_design) ---")
print(f"The design wheel load is {design_load_ton} tons, which is equivalent to {load_design_N:.2f} N.")
print(f"The tyre pressure is {tyre_pressure_kPa} kPa, which is {pressure_design_MPa:.2f} MPa.")
print("Using the formula a_design = sqrt((Force / Pressure) / π):")
print(f"a_design = sqrt(({load_design_N:.2f} N / {pressure_design_MPa:.2f} MPa) / π)")
print(f"The radius of the design wheel load is a_design = {a_design_mm:.2f} mm.\n")

# Step 3 & 4: Determine h/a ratio
# For the trial section, we find the h/a ratio
h_a_ratio_trial = h_trial_mm / a_plate_mm

# Now, we find the required F2 factor for the final design
# Δ = (p * a / E2) * F2. We can use p = Force / (π*a^2)
# Δ = ( (Force / (π*a^2)) * a / E2 ) * F2 = (Force / (π * a * E2)) * F2
# Or more simply, using the pressure: Δ = (p * a / E2) * F2
# F2_design = (Δ_design * E2) / (p_design * a_design)
F2_design = (max_deflection_mm * E2_MPa) / (pressure_design_MPa * a_design_mm)

# To verify, we can also calculate F2 for the trial pavement
# F2_trial = (Δ_trial * E2) / (p_plate * a_plate) where p_plate = P/(π*a^2)
p_plate_MPa = load_plate_N / (math.pi * a_plate_mm**2)
F2_trial = (deflection_trial_mm * E2_MPa) / (p_plate_MPa * a_plate_mm)

print("--- Step 3 & 4: Determine Thickness-to-Radius Ratio (h/a) ---")
print(f"The h/a ratio for the trial pavement is {h_trial_mm:.1f} mm / {a_plate_mm:.1f} mm = {h_a_ratio_trial:.4f}.")
print(f"The required deflection factor for the design (F2_design) is {F2_design:.4f}.")
print(f"The deflection factor from the trial section (F2_trial) is {F2_trial:.4f}.")
print("Since F2_design ≈ F2_trial, we can assume the (h/a) ratio is constant.\n")

# Step 5: Calculate Final Pavement Thickness
# h_design / a_design = h_trial / a_trial
# h_design = (h_trial / a_trial) * a_design
h_design_mm = h_a_ratio_trial * a_design_mm

print("--- Step 5: Calculate Final Pavement Thickness (h_design) ---")
print("The final pavement thickness is determined by the equation:")
print("h_design = (h_trial / a_plate) * a_design")
print("Plugging in the numbers:")
print(f"h_design = ({h_trial_mm:.1f} mm / {a_plate_mm:.1f} mm) * {a_design_mm:.2f} mm")
print(f"h_design = {h_design_mm:.2f} mm")

final_answer = h_design_mm