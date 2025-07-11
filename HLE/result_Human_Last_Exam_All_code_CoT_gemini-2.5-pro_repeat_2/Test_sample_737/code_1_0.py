import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory
    based on plate bearing test data and design load parameters.
    """
    # --- Step 1: Define given parameters and convert units ---
    # Plate bearing test on subgrade
    delta_subgrade_um = 2460.0  # μm
    load_test_kN = 30.0        # kN
    plate_diameter_mm = 305.0  # mm

    # Plate bearing test on trial pavement section
    delta_pavement_um = 1080.0 # μm
    h_test_mm = 300.0          # mm

    # Design parameters
    design_wheel_load_ton = 1.80  # metric ton
    design_pressure_kPa = 600.0   # kN/m^2
    design_deflection_limit_mm = 1.00 # mm
    
    # Constants and conversions
    g = 9.81  # m/s^2, acceleration due to gravity
    
    # Convert to consistent units (N, mm, MPa)
    delta_subgrade = delta_subgrade_um / 1000.0  # mm
    P_test = load_test_kN * 1000.0              # N
    a_plate = plate_diameter_mm / 2.0           # mm
    
    delta_pavement = delta_pavement_um / 1000.0 # mm
    h_test = h_test_mm                          # mm
    
    # Design load conversion (assuming 1 ton = 1000 kg)
    P_design = design_wheel_load_ton * 1000.0 * g # N
    # Design pressure conversion (1 kPa = 10^-3 N/mm^2)
    p_design = design_pressure_kPa / 1000.0      # N/mm^2 (MPa)

    print("--- Step 1: Input Parameters (in consistent units N, mm) ---")
    print(f"Test Load (P_test): {P_test:.2f} N")
    print(f"Test Plate Radius (a_plate): {a_plate:.2f} mm")
    print(f"Subgrade Deflection (Δ_s): {delta_subgrade:.3f} mm")
    print(f"Pavement Deflection (Δ_p): {delta_pavement:.3f} mm")
    print(f"Test Pavement Thickness (h_test): {h_test:.2f} mm")
    print(f"Design Load (P_design): {P_design:.2f} N")
    print(f"Design Tyre Pressure (p_design): {p_design:.3f} N/mm^2")
    print(f"Design Deflection Limit (Δ_design): {design_deflection_limit_mm:.2f} mm")
    print("-" * 50)

    # --- Step 2: Calculate Radius of Design Load ---
    # The tire contact area is assumed to be a circle.
    # p_design = P_design / (pi * a_design^2) => a_design = sqrt(P_design / (pi * p_design))
    a_design = math.sqrt(P_design / (math.pi * p_design))
    print("--- Step 2: Calculate Radius of Design Wheel Load ---")
    print(f"a_design = sqrt(P_design / (π * p_design))")
    print(f"a_design = sqrt({P_design:.2f} / (π * {p_design:.3f})) = {a_design:.2f} mm")
    print("-" * 50)

    # --- Step 3: Calculate Subgrade Modulus (Es) ---
    # Using the formula for a flexible plate on an elastic half-space with μ=0.5: Δ = 1.5 * p * a / Es
    # The pressure under the test plate is p_test = P_test / (pi * a_plate^2)
    p_test = P_test / (math.pi * a_plate**2)
    E_s = 1.5 * p_test * a_plate / delta_subgrade
    print("--- Step 3: Calculate Subgrade Modulus (E_s) ---")
    print(f"Using formula: Δ_s = 1.5 * p_test * a_plate / E_s")
    print(f"E_s = 1.5 * p_test * a_plate / Δ_s")
    print(f"where p_test = {P_test:.2f} / (π * {a_plate:.2f}^2) = {p_test:.4f} N/mm^2")
    print(f"E_s = (1.5 * {p_test:.4f} * {a_plate:.2f}) / {delta_subgrade:.3f} = {E_s:.2f} MPa")
    print("-" * 50)

    # --- Step 4: Calculate Deflection Factor (F2) for the Test Section ---
    # F2 is the ratio of deflection on the two-layer system to the deflection on the subgrade alone.
    F2_test = delta_pavement / delta_subgrade
    print("--- Step 4: Calculate Deflection Factor (F2) for Test Section ---")
    print(f"F2_test = Δ_p / Δ_s = {delta_pavement:.3f} / {delta_subgrade:.3f} = {F2_test:.4f}")
    print("-" * 50)
    
    # --- Step 5: Calculate Required Deflection Factor (F2) for the Design ---
    # First, find the deflection that the design load would cause on the subgrade alone.
    delta_subgrade_design = 1.5 * p_design * a_design / E_s
    # Then, find the required F2 factor to meet the design limit.
    F2_design = design_deflection_limit_mm / delta_subgrade_design
    print("--- Step 5: Calculate Required F2 for Design ---")
    print(f"Theoretical deflection on subgrade from design load: Δ_s_design = 1.5 * p_design * a_design / E_s")
    print(f"Δ_s_design = (1.5 * {p_design:.3f} * {a_design:.2f}) / {E_s:.2f} = {delta_subgrade_design:.3f} mm")
    print(f"Required F2_design = Δ_design_limit / Δ_s_design = {design_deflection_limit_mm:.2f} / {delta_subgrade_design:.3f} = {F2_design:.4f}")
    print("-" * 50)

    # --- Step 6: Determine Required Pavement Thickness (h_design) ---
    # In Burmister's theory, F2 is a function of h/a and Ep/Es.
    # Since F2_test is almost identical to F2_design, and Ep/Es is constant, we can assume h/a is also constant.
    h_design = h_test * (a_design / a_plate)
    print("--- Step 6: Final Calculation for Design Thickness (h_design) ---")
    print("Since F2_test is almost identical to F2_design (≈0.439), we assume the ratio (h/a) is the same for both cases.")
    print("h_design / a_design = h_test / a_plate")
    print("Therefore, h_design = h_test * (a_design / a_plate)")
    print(f"h_design = {h_test:.1f} * ({a_design:.2f} / {a_plate:.1f})")
    print(f"The required pavement thickness is {h_design:.2f} mm.")

    return h_design

final_thickness = solve_pavement_thickness()
print(f"\n<<<>>>")
print(f"h_design = {300.0:.1f} * ({96.78:.2f} / {152.5:.1f}) = {final_thickness:.2f} mm")
print(f"<<<{final_thickness:.2f}>>>")