import math

# --- 1. Define constants and given values in a consistent unit system (N, mm, MPa) ---

# Plate bearing test parameters
load_plate_kN = 30.0  # kN
diameter_plate_mm = 305.0  # mm
deflection_subgrade_um = 2460.0  # μm
deflection_trial_um = 1080.0  # μm
thickness_trial_mm = 300.0  # mm

# Design wheel load parameters
load_wheel_ton = 1.80  # metric ton
pressure_wheel_kPa = 600.0  # kN/m^2 or kPa
deflection_design_limit_mm = 1.00 # mm

# Conversions
g = 9.81 # m/s^2
load_plate_N = load_plate_kN * 1000.0
load_wheel_N = load_wheel_ton * 1000.0 * g
pressure_wheel_MPa = pressure_wheel_kPa * 0.001 # 1 kPa = 0.001 MPa (N/mm^2)
deflection_subgrade_mm = deflection_subgrade_um / 1000.0
deflection_trial_mm = deflection_trial_um / 1000.0
radius_plate_mm = diameter_plate_mm / 2.0

print("--- Input Parameters (Converted to N, mm) ---")
print(f"Plate Load (P_plate): {load_plate_N:.2f} N")
print(f"Plate Radius (a_plate): {radius_plate_mm:.2f} mm")
print(f"Subgrade Deflection (Δ_s): {deflection_subgrade_mm:.3f} mm")
print(f"Trial Pavement Deflection (Δ_p): {deflection_trial_mm:.3f} mm")
print(f"Trial Pavement Thickness (h_trial): {thickness_trial_mm:.2f} mm")
print(f"Design Wheel Load (P_design): {load_wheel_N:.2f} N")
print(f"Design Tyre Pressure (p_design): {pressure_wheel_MPa:.3f} N/mm^2")
print(f"Design Deflection Limit (Δ_design): {deflection_design_limit_mm:.2f} mm\n")

# --- 2. Calculate Design Load Radius ---
# Area = Force / Pressure
area_design_mm2 = load_wheel_N / pressure_wheel_MPa
# Area = pi * r^2  => r = sqrt(Area / pi)
radius_design_mm = math.sqrt(area_design_mm2 / math.pi)
print("--- Step 1: Calculate Design Load Radius ---")
print(f"Design load contact radius (a_design): {radius_design_mm:.2f} mm\n")

# --- 3. Determine Subgrade Modulus (E₂) ---
# Formula for deflection on elastic half-space: Δ = (1.5 * P) / (π * a * E)
# Rearranging for E: E = (1.5 * P) / (π * a * Δ)
E2_subgrade_MPa = (1.5 * load_plate_N) / (math.pi * radius_plate_mm * deflection_subgrade_mm)
print("--- Step 2: Calculate Subgrade Modulus (E₂) ---")
print(f"Subgrade Modulus (E₂): {E2_subgrade_MPa:.2f} MPa\n")

# --- 4. Analyze the Two-Layer System (Trial) ---
# F₂ is the ratio of deflection with the pavement layer to deflection without it.
# F₂ = Δ_trial / Δ_subgrade
F2_trial = deflection_trial_mm / deflection_subgrade_mm
# Calculate the thickness-to-radius ratio for the trial
h_a_ratio_trial = thickness_trial_mm / radius_plate_mm
print("--- Step 3: Analyze Trial Section ---")
print(f"Deflection Factor for trial (F₂_trial): {F2_trial:.4f}")
print(f"Thickness/Radius ratio for trial (h/a)_trial: {h_a_ratio_trial:.4f}\n")

# --- 5. Determine Required Pavement Thickness (h_design) ---
# Calculate the required F₂ to meet the design deflection limit
# Δ_design = (1.5 * P_design * F₂_design) / (π * a_design * E₂)
# Rearranging for F₂_design:
F2_design_req = (deflection_design_limit_mm * math.pi * radius_design_mm * E2_subgrade_MPa) / (1.5 * load_wheel_N)
print("--- Step 4: Determine Required Design Parameters ---")
print(f"Required Deflection Factor (F₂_design): {F2_design_req:.4f}\n")

# Since F₂_design is very close to F₂_trial, we can assume (h/a)_design ≈ (h/a)_trial
# (h_design / a_design) = (h_trial / a_plate)
# Rearranging for h_design:
h_design_mm = h_a_ratio_trial * radius_design_mm
print("--- Step 5: Calculate Final Pavement Thickness ---")
print("Since F₂_design is nearly identical to F₂_trial, we assume their (h/a) ratios are equal.")
print("The final calculation is based on the proportion: h_design / a_design = h_trial / a_plate")
print("\nFinal Equation:")
print(f"Required Pavement Thickness = ({thickness_trial_mm:.1f} mm / {radius_plate_mm:.1f} mm) * {radius_design_mm:.2f} mm")
print(f"Result: {h_design_mm:.1f} mm")

# --- Final Answer ---
print(f"\n<<<The required pavement thickness is {h_design_mm:.1f} mm>>>")