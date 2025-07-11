import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Step 0: Initial Given Data & Conversions ---
    # Plate Bearing Test on Subgrade
    delta_s_um = 2460
    P_plate_kN = 30
    d_plate_mm = 305

    # Plate Bearing Test on Trial Section
    delta_t_um = 1080
    h_trial_mm = 300

    # Design Parameters
    W_design_ton = 1.80
    q_design_kpa = 600
    delta_allowable_mm = 1.00
    
    # Constants & Unit Conversions
    g = 9.81  # m/s^2
    delta_s_mm = delta_s_um / 1000
    delta_t_mm = delta_t_um / 1000
    P_plate_N = P_plate_kN * 1000
    a_plate_mm = d_plate_mm / 2
    P_design_N = W_design_ton * 1000 * g
    q_design_mpa = q_design_kpa / 1000  # 1 MPa = 1 N/mm^2

    print("--- Problem Setup ---")
    print(f"Allowable deflection: {delta_allowable_mm} mm")
    print(f"Plate test load: {P_plate_kN} kN, Plate diameter: {d_plate_mm} mm")
    print(f"Design wheel load: {W_design_ton} ton, Tyre pressure: {q_design_kpa} kPa\n")

    # --- Step 1: Determine Subgrade Modulus (Es) ---
    # Using the formula for a rigid plate: Δ = 1.18 * P / (π * a * E)
    E_s_mpa = (1.18 * P_plate_N) / (math.pi * a_plate_mm * delta_s_mm)
    print("--- Step 1: Calculate Subgrade Modulus (Es) ---")
    print(f"Using subgrade deflection Δ_s = {delta_s_mm:.2f} mm")
    print(f"Equation: Es = (1.18 * P) / (π * a * Δ_s)")
    print(f"Es = (1.18 * {P_plate_N:.0f}) / (π * {a_plate_mm:.1f} * {delta_s_mm:.2f})")
    print(f"Result: Es = {E_s_mpa:.2f} MPa\n")

    # --- Step 2: Determine Two-Layer Deflection Factor (F2) from Trial ---
    # For a two-layer system, Δ_t = Δ_s * F2
    F2_test = delta_t_mm / delta_s_mm
    h_a_test_ratio = h_trial_mm / a_plate_mm
    print("--- Step 2: Calculate F2 from Trial Section ---")
    print(f"Using pavement deflection Δ_t = {delta_t_mm:.2f} mm")
    print(f"Equation: F2 = Δ_t / Δ_s")
    print(f"F2 = {delta_t_mm:.2f} / {delta_s_mm:.2f} = {F2_test:.4f}")
    print(f"Trial section's h/a ratio = {h_trial_mm:.0f} mm / {a_plate_mm:.1f} mm = {h_a_test_ratio:.4f}\n")

    # --- Step 3: Find Modular Ratio (Ep/Es) by Interpolation ---
    # Using Burmister chart data for h/a ≈ 1.967 (approximated as 2.0).
    # At h/a = 2.0: (Ep/Es=10, F2=0.48), (Ep/Es=20, F2=0.35)
    f2_1, epes_1 = 0.48, 10
    f2_2, epes_2 = 0.35, 20
    Ep_Es_ratio = epes_1 + (F2_test - f2_1) * (epes_2 - epes_1) / (f2_2 - f2_1)
    print("--- Step 3: Find Modular Ratio (Ep/Es) ---")
    print(f"Interpolating from Burmister chart data at h/a ≈ 2.0 for F2 = {F2_test:.4f}")
    print(f"Equation: Ep/Es = {epes_1} + ({F2_test:.4f} - {f2_1}) * ({epes_2} - {epes_1}) / ({f2_2} - {f2_1})")
    print(f"Result: Ep/Es = {Ep_Es_ratio:.2f}\n")

    # --- Step 4: Calculate Design Load Radius (a_design) ---
    # a = sqrt(P / (q * π))
    contact_area_mm2 = P_design_N / q_design_mpa
    a_design_mm = math.sqrt(contact_area_mm2 / math.pi)
    print("--- Step 4: Calculate Design Load Radius (a_design) ---")
    print(f"Design Load P_design = {P_design_N:.0f} N, Tyre Pressure q_design = {q_design_mpa:.3f} N/mm^2")
    print(f"Equation: a_design = sqrt(P_design / (q_design * π))")
    print(f"a_design = sqrt({P_design_N:.0f} / ({q_design_mpa:.3f} * π))")
    print(f"Result: a_design = {a_design_mm:.2f} mm\n")

    # --- Step 5: Calculate Required Design F2 (F2_design) ---
    # For a flexible tire load: Δ = (1.5 * q * a * F2) / E_s
    F2_design = (delta_allowable_mm * E_s_mpa) / (1.5 * q_design_mpa * a_design_mm)
    print("--- Step 5: Calculate Required F2 for Design ---")
    print("Using flexible load formula (factor 1.5) for the tire.")
    print(f"Equation: F2_design = (Δ_allowable * Es) / (1.5 * q_design * a_design)")
    print(f"F2_design = ({delta_allowable_mm:.2f} * {E_s_mpa:.2f}) / (1.5 * {q_design_mpa:.3f} * {a_design_mm:.2f})")
    print(f"Result: F2_design = {F2_design:.4f}\n")

    # --- Step 6: Find Design h/a Ratio by Interpolation ---
    # Interpolate to find h/a for F2_design and Ep_Es_ratio.
    # Chart data: Curve Ep/Es=10: (h/a=2, F2=0.48), (h/a=3, F2=0.35)
    #             Curve Ep/Es=20: (h/a=2, F2=0.35), (h/a=3, F2=0.25)
    h_a_on_EpEs10 = 2.0 + (F2_design - 0.48) * (3.0 - 2.0) / (0.35 - 0.48)
    h_a_on_EpEs20 = 2.0 + (F2_design - 0.35) * (3.0 - 2.0) / (0.25 - 0.35)
    h_a_design_ratio = h_a_on_EpEs10 + (Ep_Es_ratio - 10) * (h_a_on_EpEs20 - h_a_on_EpEs10) / (20 - 10)
    print("--- Step 6: Find Design h/a Ratio ---")
    print(f"Interpolating on Burmister charts for F2 = {F2_design:.4f} and Ep/Es = {Ep_Es_ratio:.2f}")
    print(f"Intermediate h/a on Ep/Es=10 curve: {h_a_on_EpEs10:.4f}")
    print(f"Intermediate h/a on Ep/Es=20 curve: {h_a_on_EpEs20:.4f}")
    print(f"Final interpolation equation for h/a_design:")
    print(f"h/a = {h_a_on_EpEs10:.4f} + ({Ep_Es_ratio:.2f} - 10) * ({h_a_on_EpEs20:.4f} - {h_a_on_EpEs10:.4f}) / (20 - 10)")
    print(f"Result: h/a_design = {h_a_design_ratio:.4f}\n")

    # --- Step 7: Calculate Final Pavement Thickness ---
    h_design_mm = h_a_design_ratio * a_design_mm
    print("--- Step 7: Calculate Final Pavement Thickness (h_design) ---")
    print(f"Equation: h_design = h/a_design * a_design")
    print(f"h_design = {h_a_design_ratio:.4f} * {a_design_mm:.2f}")
    print(f"Result: h_design = {h_design_mm:.2f} mm")

    return h_design_mm

if __name__ == '__main__':
    final_thickness = solve_pavement_thickness()
    print("\n" + "="*40)
    print(f"FINAL ANSWER:")
    print(f"The required pavement thickness is {final_thickness:.1f} mm.")
    print("="*40)
    # The final answer format for the platform
    # <<<264.2>>>