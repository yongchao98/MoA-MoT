import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Step 0: Given Data and Unit Conversions ---
    # Convert all inputs to consistent units: Newtons (N) and millimeters (mm)

    # Test 1 (Subgrade)
    w_s_microm = 2460      # Deflection (μm)
    P_s_kN = 30            # Load (kN)
    d_s_mm = 305           # Plate diameter (mm)
    w_s_mm = w_s_microm / 1000.0
    P_s_N = P_s_kN * 1000.0
    a_s_mm = d_s_mm / 2.0

    # Test 2 (Trial Pavement)
    w_p_microm = 1080      # Deflection (μm)
    h_test_mm = 300        # Pavement thickness (mm)
    w_p_mm = w_p_microm / 1000.0

    # Design Parameters
    w_design_max_mm = 1.00 # Max allowable deflection (mm)
    W_design_ton = 1.80    # Design wheel load (tons)
    p_design_kPa = 600     # Tyre pressure (kN/m^2)
    g = 9.81               # Gravitational acceleration (m/s^2)
    P_design_N = W_design_ton * 1000.0 * g
    p_design_N_mm2 = p_design_kPa / 1000.0 # 1 kPa = 0.001 N/mm^2

    # --- Step 1: Calculate Subgrade Modulus (E₂) ---
    # Pressure from plate bearing test: p = P / (pi * a^2)
    p_s_N_mm2 = P_s_N / (math.pi * a_s_mm**2)
    # Boussinesq's equation for deflection (μ = 0.5): w = 1.5 * p * a / E
    # Rearranged for E₂: E₂ = 1.5 * p_s * a_s / w_s
    E2_N_mm2 = (1.5 * p_s_N_mm2 * a_s_mm) / w_s_mm

    # --- Step 2: Analyze Trial Section ---
    # Deflection factor F₂ from trial: w_p = w_s * F₂
    F2_test = w_p_mm / w_s_mm
    # Thickness-to-radius ratio from trial: h/a
    h_a_ratio_test = h_test_mm / a_s_mm

    # --- Step 3: Determine Design Requirements ---
    # Radius of design load: P = p * A => A = P/p => pi*a^2 = P/p => a = sqrt(P/(pi*p))
    a_design_mm = math.sqrt(P_design_N / (math.pi * p_design_N_mm2))
    # Theoretical deflection if design load was on subgrade:
    w_subgrade_design_mm = (1.5 * p_design_N_mm2 * a_design_mm) / E2_N_mm2
    # Required deflection factor for design: F₂_design = w_design_max / w_subgrade_design
    F2_design = w_design_max_mm / w_subgrade_design_mm

    # --- Step 4: Calculate Required Pavement Thickness (h_design) ---
    # Since F₂_design ≈ F₂_test, we can assume (h/a)_design ≈ (h/a)_test.
    # h_design / a_design = h_test / a_s
    h_design_mm = h_a_ratio_test * a_design_mm

    # --- Final Output ---
    print("--- Analysis of the Pavement Design Problem ---\n")
    print(f"1. The modulus of the subgrade (E₂) is calculated from the first test to be {E2_N_mm2:.2f} MPa.\n")
    print(f"2. The trial section had a thickness-to-radius ratio (h_test/a_s) of {h_a_ratio_test:.4f} and produced a deflection factor (F₂) of {F2_test:.4f}.\n")
    print(f"3. The design load has an effective radius (a_design) of {a_design_mm:.2f} mm.")
    print(f"4. To limit deflection to {w_design_max_mm:.2f} mm, the required deflection factor (F₂) is calculated to be {F2_design:.4f}.\n")
    print("5. Since the required F₂ is very close to the F₂ from the trial, we assume the thickness-to-radius ratio must also be the same.")
    
    print("\n--- Final Calculation Equation ---")
    print("Required Thickness (h_design) = (h_test / a_s) * a_design")
    print("\nSubstituting the calculated values:")
    print(f"h_design = ({h_test_mm:.2f} mm / {a_s_mm:.2f} mm) * {a_design_mm:.2f} mm")
    print(f"h_design = {h_a_ratio_test:.4f} * {a_design_mm:.2f} mm")
    
    print("\n--- Result ---")
    print(f"The required pavement thickness is: {h_design_mm:.2f} mm")

    return h_design_mm

if __name__ == '__main__':
    final_thickness = solve_pavement_thickness()
    print(f"\n<<< {final_thickness:.2f} >>>")
