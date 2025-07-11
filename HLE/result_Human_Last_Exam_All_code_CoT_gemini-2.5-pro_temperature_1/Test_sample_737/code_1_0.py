import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    All calculations will be done in Newtons (N) and millimeters (mm).
    """
    print("### Pavement Thickness Calculation using Burmister's Two-Layer Theory ###\n")

    # --- Given Data ---
    # Plate Bearing Test 1 (on subgrade)
    P_plate = 30 * 1000  # Load in N
    d_plate = 305  # Plate diameter in mm
    delta_s = 2460 / 1000  # Deflection on subgrade in mm

    # Plate Bearing Test 2 (on trial pavement)
    h_trial = 300  # Trial pavement thickness in mm
    delta_t = 1080 / 1000  # Deflection on trial section in mm

    # Design Load
    W_design_kg = 1.80 * 1000  # Design wheel load in kg
    P_design = W_design_kg * 9.81 # Design wheel load in N
    p_design_kPa = 600  # Tyre pressure in kN/m^2 (kPa)
    p_design_MPa = p_design_kPa / 1000 # Tyre pressure in N/mm^2 (MPa)

    # Design Constraint
    delta_design = 1.00  # Maximum allowable deflection in mm
    mu = 0.5 # Poisson's ratio for all materials

    # --- Step 1: Determine Subgrade Modulus (E_s) ---
    print("--- Step 1: Calculating Subgrade Modulus (E_s) ---")
    a_plate = d_plate / 2
    # For a rigid plate with mu=0.5, the formula for deflection (delta) is:
    # delta = 1.178 * P / (E * a). Rearranging for E: E = 1.178 * P / (delta * a)
    # A slightly more general formula is delta = (pi/2)*(1-mu^2)*P / (E*a), which also gives E = 1.178 * P / (delta*a) for mu=0.5
    E_s = (1.178 * P_plate) / (delta_s * a_plate)
    print(f"Plate Radius (a_plate): {a_plate:.1f} mm")
    print(f"Subgrade Deflection (Δ_s): {delta_s:.3f} mm")
    print(f"Formula: E_s = (1.178 * P_plate) / (Δ_s * a_plate)")
    print(f"E_s = (1.178 * {P_plate} N) / ({delta_s:.3f} mm * {a_plate:.1f} mm)")
    print(f"Result: Subgrade Modulus (E_s) = {E_s:.2f} N/mm^2 (MPa)\n")

    # --- Step 2: Determine Pavement Modulus (E_p) ---
    print("--- Step 2: Calculating Pavement Modulus (E_p) ---")
    # Deflection factor F = delta_two_layer / delta_subgrade_only
    F_trial = delta_t / delta_s
    h_over_a_trial = h_trial / a_plate
    print(f"Deflection Factor (F_trial) = Δ_t / Δ_s = {delta_t:.3f} mm / {delta_s:.3f} mm = {F_trial:.4f}")
    print(f"Thickness/Radius Ratio (h/a_plate) = {h_trial} mm / {a_plate:.1f} mm = {h_over_a_trial:.4f}")

    # Interpolate from Burmister chart data to find Ep/Es
    # Chart data points for F at h/a ≈ 2.0:
    # (Ep/Es, F): (5, 0.48), (10, 0.35)
    # Linear interpolation: y = mx + c where y=F, x=Ep/Es
    m = (0.35 - 0.48) / (10 - 5)
    c = 0.48 - m * 5
    Ep_over_Es = (F_trial - c) / m
    print(f"Using interpolation based on Burmister's chart data for h/a ≈ 2...")
    print(f"Result: Modulus Ratio (E_p/E_s) ≈ {Ep_over_Es:.2f}")

    E_p = Ep_over_Es * E_s
    print(f"Pavement Modulus (E_p) = E_s * (E_p/E_s) = {E_s:.2f} MPa * {Ep_over_Es:.2f} = {E_p:.2f} MPa\n")

    # --- Step 3: Determine Required Pavement Thickness (h_design) ---
    print("--- Step 3: Calculating Required Pavement Thickness (h_design) ---")
    # Calculate radius of design wheel load (flexible load)
    load_area_design = P_design / p_design_MPa
    a_design = math.sqrt(load_area_design / math.pi)
    print(f"Design Load (P_design): {P_design:.2f} N")
    print(f"Design Tyre Pressure (p_design): {p_design_MPa:.3f} MPa")
    print(f"Design Load Radius (a_design) = sqrt(P_design / (p_design * pi)) = {a_design:.1f} mm")

    # Calculate the deflection on subgrade alone under the design load
    # For a flexible load (tyre) with mu=0.5, delta = 1.5 * p * a / E
    delta_subgrade_only_design = (1.5 * p_design_MPa * a_design) / E_s
    print(f"Hypothetical deflection on subgrade (Δ_subgrade_only) = (1.5 * {p_design_MPa:.3f} MPa * {a_design:.1f} mm) / {E_s:.2f} MPa = {delta_subgrade_only_design:.3f} mm")

    # Required deflection factor F_design
    F_design = delta_design / delta_subgrade_only_design
    print(f"Required Deflection Factor (F_design) = Δ_design / Δ_subgrade_only = {delta_design:.2f} mm / {delta_subgrade_only_design:.3f} mm = {F_design:.4f}")

    # Interpolate again to find h/a ratio for F_design and the known Ep_over_Es
    print(f"Using interpolation to find required (h/a) for F={F_design:.4f} at E_p/E_s={Ep_over_Es:.2f}...")
    # Chart data for h/a vs F for Ep/Es=5 and Ep/Es=10
    # For Ep/Es = 5: (h/a, F) -> (2, 0.48), (3, 0.35)
    m1 = (0.35 - 0.48) / (3 - 2)
    c1 = 0.48 - m1 * 2
    h_a_for_Ep5 = (F_design - c1) / m1

    # For Ep/Es = 10: (h/a, F) -> (2, 0.35), (3, 0.25)
    m2 = (0.25 - 0.35) / (3 - 2)
    c2 = 0.35 - m2 * 2
    h_a_for_Ep10 = (F_design - c2) / m2

    # Interpolate between the two h/a values
    # Points are (Ep/Es, h/a) -> (5, h_a_for_Ep5), (10, h_a_for_Ep10)
    m3 = (h_a_for_Ep10 - h_a_for_Ep5) / (10 - 5)
    c3 = h_a_for_Ep5 - m3 * 5
    h_over_a_design = m3 * Ep_over_Es + c3
    print(f"Result: Required Thickness/Radius Ratio (h/a)_design = {h_over_a_design:.2f}\n")

    # --- Step 4: Final Calculation ---
    print("--- Step 4: Final Pavement Thickness Calculation ---")
    h_design = h_over_a_design * a_design
    print(f"Formula: h_design = (h/a)_design * a_design")
    print(f"h_design = {h_over_a_design:.2f} * {a_design:.1f} mm")
    print(f"Final Required Pavement Thickness = {h_design:.1f} mm")
    
    return h_design

if __name__ == '__main__':
    final_thickness = solve_pavement_thickness()
    # The final answer format required by the prompt
    print(f"\n<<<{final_thickness:.1f}>>>")
