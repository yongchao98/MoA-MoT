import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Given Data ---
    # Plate Bearing Test on Subgrade
    delta_s_um = 2460  # μm
    P_plate_kN = 30  # kN
    d_plate_mm = 305  # mm

    # Plate Bearing Test on Trial Pavement
    delta_t_um = 1080  # μm
    h_trial_mm = 300  # mm

    # Design Wheel Load
    W_design_ton = 1.80  # ton
    p_design_kPa = 600  # kN/m^2

    # Design Constraint
    delta_allowable_mm = 1.00  # mm

    # --- Unit Conversions and Initial Calculations ---
    delta_s_mm = delta_s_um / 1000
    delta_t_mm = delta_t_um / 1000
    a_plate_mm = d_plate_mm / 2
    P_plate_N = P_plate_kN * 1000
    p_design_MPa = p_design_kPa / 1000
    P_design_N = W_design_ton * 1000 * 9.81  # F = m*g, assuming g = 9.81 m/s^2

    print("--- Step 1: Analyze Plate Bearing Tests ---")
    
    # Calculate deflection factor F2 from test results
    F2_trial = delta_t_mm / delta_s_mm
    print(f"1a. Deflection with trial pavement, Δ_t = {delta_t_mm} mm")
    print(f"1b. Deflection on subgrade, Δ_s = {delta_s_mm} mm")
    print(f"1c. The deflection factor F2 is the ratio Δ_t / Δ_s.")
    print(f"    F2 = {delta_t_mm} / {delta_s_mm} = {F2_trial:.4f}")
    
    # Calculate thickness-to-radius ratio for the trial section
    h_a_ratio_trial = h_trial_mm / a_plate_mm
    print(f"\n1d. The thickness-to-radius ratio for the trial is h_trial / a_plate.")
    print(f"    h_trial / a_plate = {h_trial_mm} mm / {a_plate_mm} mm = {h_a_ratio_trial:.4f}")

    print("\n1e. From Burmister's chart, for F2 ≈ 0.439 and h/a ≈ 1.97, the modular ratio E_p/E_s is approximately 10.")
    print("    We will use this established material property for the design.")

    # Calculate subgrade modulus Es
    p_plate_MPa = P_plate_N / (math.pi * a_plate_mm**2)
    E_s_MPa = (1.5 * p_plate_MPa * a_plate_mm) / delta_s_mm
    print(f"\n1f. For later steps, we calculate the subgrade modulus, E_s.")
    print(f"    Plate pressure, p_plate = {p_plate_MPa:.4f} MPa")
    print(f"    Subgrade Modulus, E_s = (1.5 * p_plate * a_plate) / Δ_s = {E_s_MPa:.2f} MPa")
    
    print("\n--- Step 2: Determine Design Parameters ---")

    # Calculate radius of design wheel contact area
    a_design_m2 = (P_design_N / 1000) / p_design_kPa
    a_design_m = math.sqrt(a_design_m2 / math.pi)
    a_design_mm = a_design_m * 1000
    print(f"2a. Design wheel load, P_design = {P_design_N/1000:.2f} kN")
    print(f"2b. Design tyre pressure, p_design = {p_design_kPa} kPa = {p_design_MPa} MPa")
    print(f"2c. Radius of design contact area, a_design = sqrt(P_design / (π * p_design)) = {a_design_mm:.2f} mm")

    # Calculate required deflection factor F2_design
    delta_subgrade_equiv = (1.5 * p_design_MPa * a_design_mm) / E_s_MPa
    F2_design = delta_allowable_mm / delta_subgrade_equiv
    print(f"\n2d. Theoretical deflection on subgrade under design load = {delta_subgrade_equiv:.4f} mm")
    print(f"2e. Required deflection factor, F2_design = Δ_allowable / Δ_subgrade_equiv")
    print(f"    F2_design = {delta_allowable_mm} / {delta_subgrade_equiv:.4f} = {F2_design:.4f}")

    print("\n--- Step 3: Calculate Required Pavement Thickness ---")
    print("3a. We use the modular ratio (E_p/E_s ≈ 10) and the required F2_design to find the required h/a ratio from the chart.")
    print(f"3b. Since F2_design ({F2_design:.4f}) is almost identical to F2_trial ({F2_trial:.4f}), the required h/a ratio for the design must be the same as the trial ratio.")
    print(f"    Therefore, h_design / a_design = h_trial / a_plate = {h_a_ratio_trial:.4f}")

    # Calculate final design thickness
    h_design_mm = h_a_ratio_trial * a_design_mm
    print("\n3c. The required pavement thickness, h_design, can now be calculated:")
    print(f"    h_design = (h_trial / a_plate) * a_design")
    print(f"    h_design = {h_a_ratio_trial:.4f} * {a_design_mm:.2f}")
    print(f"\nFinal Answer: The required pavement thickness is {h_design_mm:.2f} mm.")
    
    # Final answer in the requested format
    return h_design_mm

# Execute the function and print the final answer in the required format
final_thickness = solve_pavement_thickness()
print(f"\n<<< {final_thickness:.2f} >>>")
