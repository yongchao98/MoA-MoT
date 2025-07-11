import math

# --- Given Data ---
# Plate Bearing Test on Subgrade
deflection_subgrade_um = 2460.0
load_test_kN = 30.0
plate_diameter_mm = 305.0

# Plate Bearing Test on Trial Pavement
deflection_pavement_um = 1080.0
pavement_thickness_trial_mm = 300.0

# Design Parameters
design_wheel_load_ton = 1.80
tyre_pressure_kPa = 600.0
max_deflection_mm = 1.00

# Material Property
poisson_ratio = 0.5 # mu
g_acceleration = 9.81 # m/s^2 to convert ton (mass) to kN (force)

# --- Step-by-step Calculation ---

print("Step 1: Convert units to a consistent system (N, mm, MPa)")
# Convert deflections from micrometers to millimeters
deflection_subgrade_mm = deflection_subgrade_um / 1000.0
deflection_pavement_mm = deflection_pavement_um / 1000.0
# Convert loads from kN to N
load_test_N = load_test_kN * 1000.0
design_load_N = design_wheel_load_ton * 1000.0 * g_acceleration
# Convert pressure from kN/m^2 (kPa) to N/mm^2 (MPa)
tyre_pressure_MPa = tyre_pressure_kPa / 1000.0
print(f"Test load: {load_test_N:.2f} N")
print(f"Design load: {design_load_N:.2f} N")
print(f"Tyre pressure: {tyre_pressure_MPa:.2f} N/mm^2 (MPa)\n")


print("Step 2: Analyze the plate bearing test on the subgrade to find its elastic modulus (E2)")
# Calculate test plate radius and area
plate_radius_mm = plate_diameter_mm / 2.0
plate_area_mm2 = math.pi * plate_radius_mm**2
# Calculate pressure under the test plate
test_pressure_MPa = load_test_N / plate_area_mm2
# The formula for deflection at the center of a flexible circular load on an elastic half-space is:
# Delta = (1.5 * p * a) / E  (for Poisson's ratio mu = 0.5)
# Rearranging for E2 (subgrade modulus):
E2_subgrade_MPa = (1.5 * test_pressure_MPa * plate_radius_mm) / deflection_subgrade_mm
print(f"Test plate radius (a_test): {plate_radius_mm:.2f} mm")
print(f"Test pressure (p_test): {test_pressure_MPa:.4f} MPa")
print(f"Calculated Subgrade Modulus (E2): {E2_subgrade_MPa:.2f} MPa\n")


print("Step 3: Analyze the design wheel load to find the design radius (a_design)")
# Calculate the area of the design wheel load
design_load_area_mm2 = design_load_N / tyre_pressure_MPa
# Calculate the equivalent radius of the design load
design_radius_mm = math.sqrt(design_load_area_mm2 / math.pi)
print(f"Design load area: {design_load_area_mm2:.2f} mm^2")
print(f"Design load radius (a_design): {design_radius_mm:.2f} mm\n")


print("Step 4: Calculate Burmister's deflection factor (F2) for both trial and design cases")
# The deflection of a two-layer system is Delta = (1.5 * p * a / E2) * F2
# So, F2 = Delta_measured / Delta_on_subgrade_alone
# For the trial case:
F2_trial = deflection_pavement_mm / deflection_subgrade_mm
print(f"F2 for trial section: {F2_trial:.4f}")

# For the design case, F2 is calculated from the design parameters:
# F2_design = (Delta_design * E2) / (1.5 * p_design * a_design)
F2_design = (max_deflection_mm * E2_subgrade_MPa) / (1.5 * tyre_pressure_MPa * design_radius_mm)
print(f"Required F2 for design: {F2_design:.4f}\n")


print("Step 5: Determine the required pavement thickness (h_design)")
print("Since F2_trial is nearly identical to F2_design, we can assume their h/a ratios are equal.")
# Calculate the h/a ratio for the trial section
h_a_ratio_trial = pavement_thickness_trial_mm / plate_radius_mm
print(f"(h/a) ratio for trial section = {pavement_thickness_trial_mm:.2f} mm / {plate_radius_mm:.2f} mm = {h_a_ratio_trial:.4f}")

# Use this ratio to find the required design thickness
# (h/a)_design = (h/a)_trial
# h_design = (h/a)_trial * a_design
required_thickness_mm = h_a_ratio_trial * design_radius_mm
print("\n--- Final Calculation ---")
print(f"Required pavement thickness h = (h/a)_ratio * a_design")
print(f"h = {h_a_ratio_trial:.4f} * {design_radius_mm:.2f} mm")
print(f"h = {required_thickness_mm:.2f} mm")

<<<190.49>>>