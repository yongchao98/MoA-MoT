import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    print("### Pavement Thickness Calculation using Burmister's Two-Layer Theory ###\n")

    # --- Inputs and Constants ---
    # Test 1: On subgrade
    P_plate = 30e3  # N (30 kN)
    D_plate = 305.0  # mm
    delta_z1 = 2460.0 / 1000  # mm (converted from 2460 μm)

    # Test 2: On trial pavement section
    h_trial = 300.0  # mm
    delta_z2 = 1080.0 / 1000  # mm (converted from 1080 μm)

    # Design parameters
    design_deflection_limit = 1.00  # mm
    design_load_mass = 1.80 * 1000  # kg (1.80 ton)
    g = 9.81 # m/s^2
    P_design = design_load_mass * g # N
    p_design_kPa = 600.0  # kN/m^2
    p_design_MPa = p_design_kPa / 1000 # N/mm^2 (MPa)

    # Poisson's ratio
    mu = 0.5

    print("--- Given Information ---")
    print(f"Plate Load, P: {P_plate/1000} kN")
    print(f"Plate Diameter, D: {D_plate} mm")
    print(f"Deflection on Subgrade, Δ_z1: {delta_z1} mm")
    print(f"Trial Pavement Thickness, h_trial: {h_trial} mm")
    print(f"Deflection on Trial Section, Δ_z2: {delta_z2} mm")
    print(f"Design Load, P_design: {P_design/1000:.2f} kN")
    print(f"Tyre Pressure, p_design: {p_design_MPa} MPa")
    print(f"Max Allowable Deflection: {design_deflection_limit} mm")
    print("-" * 25 + "\n")

    # --- Step 1: Determine Subgrade Modulus (E2) ---
    print("--- Step 1: Calculating Subgrade Modulus (E₂) ---")
    a_plate = D_plate / 2
    p_plate = P_plate / (math.pi * a_plate**2)
    # Using the single-layer formula: delta_z1 = 1.5 * p * a / E2
    E2 = (1.5 * p_plate * a_plate) / delta_z1
    print(f"Plate contact radius, a_plate: {a_plate:.2f} mm")
    print(f"Plate contact pressure, p_plate: {p_plate:.4f} MPa")
    print(f"Subgrade Modulus, E₂: {E2:.2f} MPa")
    print("-" * 25 + "\n")

    # --- Step 2: Determine Pavement Modulus (E1) ---
    print("--- Step 2: Calculating Pavement Modulus (E₁) ---")
    # Using the two-layer formula: delta_z2 = (1.5 * p * a / E2) * F2
    # First, find the F2 factor from the trial section test data
    F2_trial = (delta_z2 * E2) / (1.5 * p_plate * a_plate)
    h_div_a_trial = h_trial / a_plate
    print(f"From trial section data, calculated F₂ is: {F2_trial:.4f}")
    print(f"The h/a ratio for the trial is: {h_div_a_trial:.4f}")

    # Burmister F2 chart data for mu=0.5 (as read from standard charts)
    # Format: {E1/E2_ratio: {h/a_ratio: F2_value}}
    f2_chart_data = {
        5: {1.0: 0.66, 2.0: 0.49},
        10: {1.0: 0.55, 2.0: 0.38}
    }

    # Interpolate to find E1/E2
    # We need to find the modulus ratio E1/E2 that corresponds to F2=0.4378 at h/a=1.9672
    # First, find F2 values at h/a=1.9672 for E1/E2=5 and E1/E2=10
    h_a_pts = [1.0, 2.0]
    # For E1/E2 = 5
    f2_pts_for_E5 = [f2_chart_data[5][1.0], f2_chart_data[5][2.0]]
    f2_at_trial_ha_for_E5 = f2_pts_for_E5[0] + (h_div_a_trial - h_a_pts[0]) * (f2_pts_for_E5[1] - f2_pts_for_E5[0]) / (h_a_pts[1] - h_a_pts[0])
    # For E1/E2 = 10
    f2_pts_for_E10 = [f2_chart_data[10][1.0], f2_chart_data[10][2.0]]
    f2_at_trial_ha_for_E10 = f2_pts_for_E10[0] + (h_div_a_trial - h_a_pts[0]) * (f2_pts_for_E10[1] - f2_pts_for_E10[0]) / (h_a_pts[1] - h_a_pts[0])

    # Now interpolate between the E1/E2 values to find the ratio for our F2_trial
    E_ratio_pts = [5, 10]
    f2_pts_for_E_interp = [f2_at_trial_ha_for_E5, f2_at_trial_ha_for_E10]
    E1_div_E2 = E_ratio_pts[0] + (F2_trial - f2_pts_for_E_interp[0]) * (E_ratio_pts[1] - E_ratio_pts[0]) / (f2_pts_for_E_interp[1] - f2_pts_for_E_interp[0])
    E1 = E1_div_E2 * E2
    print(f"By interpolating chart data, the modulus ratio E₁/E₂ is found to be: {E1_div_E2:.2f}")
    print(f"Pavement Modulus, E₁: {E1:.2f} MPa")
    print("-" * 25 + "\n")

    # --- Step 3: Analyze Design Load ---
    print("--- Step 3: Calculating Design Load Parameters ---")
    # P_design = p_design * pi * a_design^2 => a_design = sqrt(P_design / (pi*p_design))
    a_design = math.sqrt(P_design / (math.pi * p_design_MPa))
    print(f"Design contact radius, a_design: {a_design:.2f} mm")
    print("-" * 25 + "\n")

    # --- Step 4: Determine Required Pavement Thickness (h_design) ---
    print("--- Step 4: Determining Required Pavement Thickness (h_design) ---")
    # Using delta = (1.5 * p * a / E2) * F2, solve for the required F2
    F2_design = (design_deflection_limit * E2) / (1.5 * p_design_MPa * a_design)
    print(f"To limit deflection to {design_deflection_limit:.2f} mm, the required F₂ is: {F2_design:.4f}")

    # Now find h/a for this F2_design, using our known E1/E2 ratio. Inverse interpolation.
    # First, find F2 values at h/a=1 and h/a=2 for our specific E1/E2 ratio
    # For h/a = 1.0
    f2_pts_for_ha1 = [f2_chart_data[5][1.0], f2_chart_data[10][1.0]]
    f2_at_Eratio_for_ha1 = f2_pts_for_ha1[0] + (E1_div_E2 - E_ratio_pts[0]) * (f2_pts_for_ha1[1] - f2_pts_for_ha1[0]) / (E_ratio_pts[1] - E_ratio_pts[0])
    # For h/a = 2.0
    f2_pts_for_ha2 = [f2_chart_data[5][2.0], f2_chart_data[10][2.0]]
    f2_at_Eratio_for_ha2 = f2_pts_for_ha2[0] + (E1_div_E2 - E_ratio_pts[0]) * (f2_pts_for_ha2[1] - f2_pts_for_ha2[0]) / (E_ratio_pts[1] - E_ratio_pts[0])
    
    # We now have two points on our specific curve: (h/a, F2)
    # Point 1: (1.0, f2_at_Eratio_for_ha1)
    # Point 2: (2.0, f2_at_Eratio_for_ha2)
    # Interpolate to find the h/a_design that gives F2_design
    h_a_design_pts = [1.0, 2.0]
    f2_design_pts = [f2_at_Eratio_for_ha1, f2_at_Eratio_for_ha2]
    h_div_a_design = h_a_design_pts[0] + (F2_design - f2_design_pts[0]) * (h_a_design_pts[1] - h_a_design_pts[0]) / (f2_design_pts[1] - f2_design_pts[0])
    print(f"The required h/a ratio is found to be: {h_div_a_design:.4f}")

    # Final calculation for thickness
    h_design = h_div_a_design * a_design
    print("\n--- Final Calculation ---")
    print(f"Required Thickness h = (h/a)_design * a_design")
    print(f"Required Thickness h = {h_div_a_design:.4f} * {a_design:.2f} mm")
    print(f"Required Pavement Thickness (h): {h_design:.2f} mm")
    print("-" * 25 + "\n")
    
    return h_design

if __name__ == "__main__":
    final_thickness = solve_pavement_thickness()
    print(f"<<<{final_thickness:.2f}>>>")