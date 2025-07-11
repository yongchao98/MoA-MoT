import numpy as np

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Input Parameters ---
    # Plate bearing test on subgrade
    P_plate_kn = 30.0  # Load in kN
    d_plate_mm = 305.0 # Plate diameter in mm
    delta_subgrade_um = 2460.0 # Deflection on subgrade in micrometers

    # Plate bearing test on trial pavement
    h_trial_mm = 300.0 # Trial pavement thickness in mm
    delta_pavement_um = 1080.0 # Deflection on trial pavement in micrometers

    # Design wheel load
    W_design_ton = 1.80 # Wheel load in tons
    p_design_kpa = 600.0 # Tyre pressure in kN/m^2

    # Design constraint
    delta_target_mm = 1.00 # Maximum allowable deflection in mm
    
    # Material property
    mu = 0.5 # Poisson's ratio for all materials

    # --- Unit Conversions ---
    P_plate_n = P_plate_kn * 1000 # N
    a_plate_mm = d_plate_mm / 2.0 # mm
    delta_subgrade_mm = delta_subgrade_um / 1000.0 # mm
    delta_pavement_mm = delta_pavement_um / 1000.0 # mm
    P_design_n = W_design_ton * 1000 * 9.81 # N (assuming 1 ton = 1000 kg and g = 9.81 m/s^2)
    p_design_mpa = p_design_kpa / 1000.0 # MPa (N/mm^2)

    print("--- Step 1: Analyze Plate Bearing Tests to Find Material Properties ---")

    # Calculate subgrade modulus E2
    # Formula: delta = (1.5 * P) / (pi * a * E) => E = (1.5 * P) / (pi * a * delta)
    E2_mpa = (1.5 * P_plate_n) / (np.pi * a_plate_mm * delta_subgrade_mm)
    print(f"1a. Subgrade Modulus (E₂):")
    print(f"    E₂ = (1.5 * {P_plate_n:.0f} N) / (π * {a_plate_mm:.2f} mm * {delta_subgrade_mm:.3f} mm)")
    print(f"    E₂ = {E2_mpa:.2f} MPa\n")

    # Calculate F2 factor and h/a for the trial pavement
    F2_trial = delta_pavement_mm / delta_subgrade_mm
    h_a_ratio_trial = h_trial_mm / a_plate_mm
    print(f"1b. Trial Pavement Deflection Factor (F₂) and h/a Ratio:")
    print(f"    F₂_trial = Deflection_Pavement / Deflection_Subgrade = {delta_pavement_mm:.3f} mm / {delta_subgrade_mm:.3f} mm = {F2_trial:.4f}")
    print(f"    (h/a)_trial = {h_trial_mm:.1f} mm / {a_plate_mm:.2f} mm = {h_a_ratio_trial:.4f}\n")

    # --- Burmister Chart Data (for mu=0.5) ---
    # A digitized representation of the standard chart.
    h_a_points = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    E1_E2_ratios = np.array([1, 2, 5, 10, 20, 50, 100])
    F2_data = { # F2 values for each E1/E2 curve at the h_a_points
        1:   np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]),
        2:   np.array([0.98, 0.90, 0.85, 0.80, 0.76, 0.61, 0.42, 0.33, 0.27, 0.24]),
        5:   np.array([0.92, 0.77, 0.67, 0.60, 0.55, 0.43, 0.30, 0.23, 0.19, 0.16]),
        10:  np.array([0.88, 0.69, 0.57, 0.50, 0.45, 0.33, 0.23, 0.18, 0.15, 0.13]),
        20:  np.array([0.82, 0.62, 0.50, 0.43, 0.38, 0.27, 0.18, 0.14, 0.12, 0.10]),
        50:  np.array([0.77, 0.56, 0.44, 0.37, 0.32, 0.22, 0.15, 0.11, 0.09, 0.08]),
        100: np.array([0.74, 0.52, 0.40, 0.33, 0.29, 0.20, 0.13, 0.10, 0.08, 0.07])
    }

    # Interpolate to find E1/E2
    f2_at_ha_trial = [np.interp(h_a_ratio_trial, h_a_points, F2_data[ratio]) for ratio in E1_E2_ratios]
    
    # np.interp requires x-coordinates to be increasing. F2 decreases as E1/E2 increases.
    # So, we pass F2 values (x) and E1/E2 ratios (y) in reversed order.
    E1_E2_val = np.interp(F2_trial, f2_at_ha_trial[::-1], E1_E2_ratios[::-1])
    print(f"1c. Pavement/Subgrade Modulus Ratio (E₁/E₂) from Burmister Chart:")
    print(f"    By interpolating chart data for (h/a)={h_a_ratio_trial:.4f} and F₂={F2_trial:.4f}, we get:")
    print(f"    E₁/E₂ = {E1_E2_val:.2f}\n")
    
    print("--- Step 2: Analyze Design Wheel Load ---")
    
    # Calculate design load contact radius 'a'
    # Formula: p = P / (pi * a^2) => a = sqrt(P / (pi * p))
    a_design_mm = np.sqrt(P_design_n / (np.pi * p_design_mpa))
    print("2a. Design Load Contact Radius (a_design):")
    print(f"    a_design = √({P_design_n:.0f} N / (π * {p_design_mpa:.3f} MPa))")
    print(f"    a_design = {a_design_mm:.2f} mm\n")

    print("--- Step 3: Determine Required Pavement Thickness ---")
    
    # Calculate the required F2 factor for the design load
    # Formula: delta = (1.5 * p * a / E2) * F2 => F2 = (delta * E2) / (1.5 * p * a)
    F2_new = (delta_target_mm * E2_mpa) / (1.5 * p_design_mpa * a_design_mm)
    print("3a. Required Deflection Factor (F₂_new):")
    print(f"    F₂_new = ({delta_target_mm:.2f} mm * {E2_mpa:.2f} MPa) / (1.5 * {p_design_mpa:.3f} MPa * {a_design_mm:.2f} mm)")
    print(f"    F₂_new = {F2_new:.4f}\n")

    # Interpolate to find the required h/a ratio
    # First, create the F2 curve for our specific E1/E2 value
    e1_e2_low = E1_E2_ratios[E1_E2_ratios <= E1_E2_val][-1]
    e1_e2_high = E1_E2_ratios[E1_E2_ratios >= E1_E2_val][0]

    if e1_e2_low == e1_e2_high:
        f2_curve_new = F2_data[e1_e2_low]
    else:
        f2_curve_low = F2_data[e1_e2_low]
        f2_curve_high = F2_data[e1_e2_high]
        f2_curve_new = np.interp(E1_E2_val, [e1_e2_low, e1_e2_high], [f2_curve_low, f2_curve_high])
        
    # Interpolate on the new curve to find h/a
    # np.interp requires x-coordinates (F2) to be increasing, so we reverse lists.
    h_a_ratio_new = np.interp(F2_new, f2_curve_new[::-1], h_a_points[::-1])
    
    print("3b. Required Thickness/Radius Ratio ((h/a)_new):")
    print(f"    Using the Burmister chart for E₁/E₂ = {E1_E2_val:.2f} and F₂ = {F2_new:.4f}, we find:")
    print(f"    (h/a)_new = {h_a_ratio_new:.4f}\n")

    # Calculate final pavement thickness
    h_new_mm = h_a_ratio_new * a_design_mm
    
    print("--- Final Result: Required Pavement Thickness ---")
    print("The required pavement thickness (h) is calculated as:")
    print(f"h = (h/a)_new * a_design")
    print(f"h = {h_a_ratio_new:.4f} * {a_design_mm:.2f} mm")
    print(f"h = {h_new_mm:.2f} mm")
    print("\nTherefore, the required pavement thickness is approximately " + f"{h_new_mm:.1f} mm.")
    return h_new_mm

if __name__ == '__main__':
    final_thickness = solve_pavement_thickness()
    # The final answer format as requested by the user prompt
    print(f"\n<<< {final_thickness:.1f} >>>")
