import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Given Data ---
    # Test 1 (on subgrade)
    P_plate = 30.0  # kN
    d_plate = 305.0  # mm
    delta_s = 2460.0 / 1000  # Convert μm to mm -> 2.460 mm

    # Test 2 (on trial pavement)
    h_trial = 300.0  # mm
    delta_t = 1080.0 / 1000  # Convert μm to mm -> 1.080 mm

    # Design Load
    W_design_ton = 1.80 # ton
    P_design = W_design_ton * 9.80665 # Convert ton (force) to kN
    q_design_kpa = 600.0  # kN/m^2

    # Design Criteria
    delta_max = 1.00  # mm

    # Constants
    PI = math.pi
    a_plate = d_plate / 2.0  # mm
    q_design_mpa = q_design_kpa / 1000 # Convert kN/m^2 to kN/mm^2 (MPa)

    print("--- Step 1: Calculate Subgrade Modulus (E2) ---")
    print("The plate bearing test uses a rigid plate. The deflection on the subgrade is given by:")
    print("delta_s = (1.18 * P_plate) / (pi * a_plate * E2)")
    print("Rearranging for E2:")
    print("E2 = (1.18 * P_plate) / (pi * a_plate * delta_s)")
    
    E2 = (1.18 * P_plate) / (PI * a_plate * delta_s)
    
    print(f"E2 = (1.18 * {P_plate:.2f} kN) / ({PI:.4f} * {a_plate:.2f} mm * {delta_s:.3f} mm)")
    print(f"E2 = {E2:.4f} kN/mm^2 (or GPa)\n")

    # --- Burmister Chart Data (for interpolation) ---
    # Fw values for different h/a and E1/E2 ratios
    # Source: Standard Pavement Engineering Charts
    Fw_data = {
        # h/a ratio
        2.0: {10: 0.495, 20: 0.345},  # E1/E2 ratio: Fw value
        3.0: {10: 0.350, 20: 0.245}
    }

    print("--- Step 2: Determine Pavement Modulus (E1) ---")
    print("First, calculate the deflection factor (Fw) for the trial section.")
    print("delta_t = (1.18 * P_plate) / (pi * a_plate * E2) * Fw")
    print("Fw = (delta_t * pi * a_plate * E2) / (1.18 * P_plate)")
    
    Fw_trial = (delta_t * PI * a_plate * E2) / (1.18 * P_plate)
    
    print(f"Fw = ({delta_t:.3f} mm * {PI:.4f} * {a_plate:.2f} mm * {E2:.4f} GPa) / (1.18 * {P_plate:.2f} kN)")
    print(f"Fw = {Fw_trial:.4f}\n")

    print("Next, find the modulus ratio E1/E2.")
    print("The thickness-to-radius ratio for the trial section is h_trial / a_plate:")
    ha_ratio_trial = h_trial / a_plate
    print(f"h/a = {h_trial:.1f} mm / {a_plate:.2f} mm = {ha_ratio_trial:.4f}\n")
    
    # Interpolate to find E1/E2
    # This performs bilinear interpolation based on the chart data
    ha1, ha2 = 2.0, 3.0
    # Interpolate Fw values at h/a = ha_ratio_trial for E1/E2=10 and E1/E2=20
    fw_at_ha_for_E10 = Fw_data[ha1][10] + (ha_ratio_trial - ha1) * (Fw_data[ha2][10] - Fw_data[ha1][10]) / (ha2 - ha1)
    fw_at_ha_for_E20 = Fw_data[ha1][20] + (ha_ratio_trial - ha1) * (Fw_data[ha2][20] - Fw_data[ha1][20]) / (ha2 - ha1)
    # Interpolate E1/E2 using the calculated Fw_trial
    E1_E2_ratio = 10 + (Fw_trial - fw_at_ha_for_E10) * (20 - 10) / (fw_at_ha_for_E20 - fw_at_ha_for_E10)
    
    print("By interpolating from Burmister chart data for Fw and h/a:")
    print(f"Modulus Ratio (E1/E2) = {E1_E2_ratio:.2f}\n")

    print("Finally, calculate E1:")
    print("E1 = E1/E2 * E2")
    E1 = E1_E2_ratio * E2
    print(f"E1 = {E1_E2_ratio:.2f} * {E2:.4f} GPa = {E1:.4f} GPa\n")

    print("--- Step 3: Analyze Design Load and Required Fw ---")
    print("The design wheel load is a flexible plate. First, find its radius (a_design).")
    print("P_design = q_design * pi * a_design^2")
    print("a_design = sqrt(P_design / (pi * q_design))")
    
    a_design = math.sqrt(P_design / (PI * q_design_mpa))
    
    print(f"a_design = sqrt({P_design:.2f} kN / ({PI:.4f} * {q_design_mpa:.6f} kN/mm^2))")
    print(f"a_design = {a_design:.2f} mm\n")

    print("Next, calculate the required design deflection factor (Fw_design).")
    print("delta_max = (1.5 * q_design * a_design) / E2 * Fw_design")
    print("Fw_design = (delta_max * E2) / (1.5 * q_design * a_design)")
    
    Fw_design = (delta_max * E2) / (1.5 * q_design_mpa * a_design)
    
    print(f"Fw_design = ({delta_max:.2f} mm * {E2:.4f} GPa) / (1.5 * {q_design_mpa:.6f} kN/mm^2 * {a_design:.2f} mm)")
    print(f"Fw_design = {Fw_design:.4f}\n")

    print("--- Step 4: Determine Required Pavement Thickness (h_design) ---")
    print("Using the known E1/E2 ratio and the required Fw_design, find the required h/a ratio from the chart.")
    
    # Interpolate to find h/a
    # This is inverse bilinear interpolation
    # Interpolate Fw values at E1/E2=E1_E2_ratio for h/a=2 and h/a=3
    fw_at_Eratio_for_ha2 = Fw_data[ha1][10] + (E1_E2_ratio - 10) * (Fw_data[ha1][20] - Fw_data[ha1][10]) / (20 - 10)
    fw_at_Eratio_for_ha3 = Fw_data[ha2][10] + (E1_E2_ratio - 10) * (Fw_data[ha2][20] - Fw_data[ha2][10]) / (20 - 10)
    # Interpolate h/a using the Fw_design
    ha_ratio_design = ha1 + (Fw_design - fw_at_Eratio_for_ha2) * (ha2 - ha1) / (fw_at_Eratio_for_ha3 - fw_at_Eratio_for_ha2)
    
    print(f"Required h/a ratio = {ha_ratio_design:.4f}\n")

    print("Finally, calculate the required pavement thickness h_design.")
    print("h_design = (h/a)_design * a_design")
    h_design = ha_ratio_design * a_design
    print(f"h_design = {ha_ratio_design:.4f} * {a_design:.2f} mm")
    print(f"h_design = {h_design:.2f} mm")
    
    return h_design

if __name__ == '__main__':
    final_thickness = solve_pavement_thickness()
    # The final answer is wrapped in <<<>>>
    print(f"\n<<<The required pavement thickness is {final_thickness:.1f} mm>>>")
    # A more machine-readable format for the final answer
    print(f"<<<{final_thickness:.1f}>>>")