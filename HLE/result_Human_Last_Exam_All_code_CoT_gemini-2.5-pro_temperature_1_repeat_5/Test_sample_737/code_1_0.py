import numpy as np

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # --- Given Data ---
    # Test 1 (on subgrade)
    P_test = 30e3  # Load (N)
    d_plate = 305  # Plate diameter (mm)
    a_plate = d_plate / 2  # Plate radius (mm)
    delta_1 = 2460e-3  # Deflection on subgrade (mm)

    # Test 2 (on pavement trial section)
    h_trial = 300  # Trial pavement thickness (mm)
    delta_2 = 1080e-3  # Deflection on pavement (mm)

    # Design parameters
    W_design_ton = 1.80 # Design wheel load (ton)
    P_design = W_design_ton * 1000 * 9.81 # Design wheel load (N), g=9.81 m/s^2
    p_design_kpa = 600  # Tyre pressure (kN/m^2 or kPa)
    p_design_mpa = p_design_kpa * 1e-3 # Tyre pressure (MPa or N/mm^2)
    delta_allowable = 1.00  # Allowable deflection (mm)

    # Material properties
    mu = 0.5  # Poisson's ratio for all materials

    # Burmister's two-layer deflection factor (F2) chart data for mu=0.5
    # F2_chart[E1/E2][h/a] = F2
    F2_chart = {
        5: {1: 0.67, 2: 0.53, 4: 0.36, 8: 0.22},
        10: {1: 0.58, 2: 0.38, 4: 0.25, 8: 0.15},
    }

    print("Step 1: Calculate Subgrade Modulus (E2)")
    # Using the formula for a rigid plate on an elastic half-space:
    # Delta = (P * (1 - mu^2)) / (2 * a * E)
    # E2 = (P * (1 - mu^2)) / (2 * a * Delta)
    E2 = (P_test * (1 - mu**2)) / (2 * a_plate * delta_1)
    print(f"E2 = ({P_test:.0f} N * (1 - {mu}**2)) / (2 * {a_plate:.1f} mm * {delta_1:.3f} mm)")
    print(f"   = {E2:.2f} MPa (N/mm^2)\n")


    print("Step 2: Determine the Modular Ratio (E1/E2)")
    # For a two-layer system, Delta_2 = Delta_1_theoretical * F2
    # Where Delta_1_theoretical is the deflection on the subgrade alone.
    # We can use the measured Delta_1 as a very close approximation.
    # F2 = Delta_2 / Delta_1
    F2_test = delta_2 / delta_1
    h_a_test = h_trial / a_plate
    print(f"The deflection factor F2 from the test is calculated as:")
    print(f"F2 = Deflection on Pavement / Deflection on Subgrade = {delta_2:.3f} mm / {delta_1:.3f} mm = {F2_test:.3f}")
    print(f"The thickness/radius ratio for the test is h/a = {h_trial:.1f} mm / {a_plate:.1f} mm = {h_a_test:.3f}\n")
    
    print("Interpolating from Burmister chart data to find E1/E2 for h/a ~ 2.0 and F2 = 0.439...")
    # We will interpolate along the h/a = 2 line to find E1/E2
    ha_key = 2
    e_ratio_keys = sorted(F2_chart.keys())
    # Find E1/E2 keys that bracket our F2 value
    e1_key, e2_key = 0, 0
    for i in range(len(e_ratio_keys) - 1):
        f2_val1 = F2_chart[e_ratio_keys[i]][ha_key]
        f2_val2 = F2_chart[e_ratio_keys[i+1]][ha_key]
        if (f2_val1 >= F2_test >= f2_val2) or (f2_val1 <= F2_test <= f2_val2):
            e1_key = e_ratio_keys[i]
            e2_key = e_ratio_keys[i+1]
            break

    F2_at_e1 = F2_chart[e1_key][ha_key]
    F2_at_e2 = F2_chart[e2_key][ha_key]

    # Linear interpolation: y = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
    # Here, y is E1/E2 and x is F2
    E1_E2_ratio = e1_key + (F2_test - F2_at_e1) * (e2_key - e1_key) / (F2_at_e2 - F2_at_e1)
    E1 = E1_E2_ratio * E2
    print(f"E1/E2 = {e1_key} + ({F2_test:.3f} - {F2_at_e1}) * ({e2_key} - {e1_key}) / ({F2_at_e2} - {F2_at_e1})")
    print(f"      = {E1_E2_ratio:.2f}")
    print(f"Thus, the pavement modulus E1 = {E1_E2_ratio:.2f} * {E2:.2f} MPa = {E1:.2f} MPa\n")


    print("Step 3: Analyze the Design Load")
    # For a circular load area, P = p * A = p * pi * a^2
    # a_design = sqrt(P / (pi * p))
    a_design = np.sqrt(P_design / (np.pi * p_design_mpa))
    print(f"The design wheel load is {W_design_ton} tons = {P_design:.0f} N")
    print(f"The tyre pressure is {p_design_kpa} kPa = {p_design_mpa:.3f} MPa")
    print(f"The radius of the design load contact area, a_design:")
    print(f"a_design = sqrt({P_design:.0f} N / (pi * {p_design_mpa:.3f} MPa)) = {a_design:.2f} mm\n")


    print("Step 4: Determine Required Pavement Thickness (h)")
    # For a flexible load, Delta = (1.5 * p * a / E2) * F2
    # So, F2_design = (Delta_allowable * E2) / (1.5 * p_design * a_design)
    F2_design = (delta_allowable * E2) / (1.5 * p_design_mpa * a_design)
    print(f"First, calculate the required deflection factor F2 for the design case:")
    print(f"F2_design = ({delta_allowable:.2f} mm * {E2:.2f} MPa) / (1.5 * {p_design_mpa:.3f} MPa * {a_design:.2f} mm)")
    print(f"          = {F2_design:.3f}\n")

    print(f"Now, find the required h/a ratio for E1/E2 = {E1_E2_ratio:.2f} and F2 = {F2_design:.3f}.")
    print("This requires 2D interpolation. We will create a curve for E1/E2 = 8.03 and find h/a.")
    
    # Interpolate F2 values for our E1/E2 ratio at h/a=2 and h/a=4
    weight = (E1_E2_ratio - e1_key) / (e2_key - e1_key)
    
    F2_at_ha2 = F2_chart[e1_key][2] + weight * (F2_chart[e2_key][2] - F2_chart[e1_key][2])
    F2_at_ha4 = F2_chart[e1_key][4] + weight * (F2_chart[e2_key][4] - F2_chart[e1_key][4])
    
    # Interpolate for h/a using a log-linear relationship, which is more accurate
    # log(h/a) = log(y1) + (x - x1) * (log(y2) - log(y1)) / (x2 - x1)
    # Here, y is h/a and x is F2
    log_ha = np.log10(2) + (F2_design - F2_at_ha2) * (np.log10(4) - np.log10(2)) / (F2_at_ha4 - F2_at_ha2)
    h_a_design = 10**log_ha
    h_design = h_a_design * a_design
    
    print(f"Interpolating log(h/a) vs F2 between points (F2={F2_at_ha2:.3f}, h/a=2) and (F2={F2_at_ha4:.3f}, h/a=4):")
    print(f"log(h/a)_design = log10(2) + ({F2_design:.3f} - {F2_at_ha2:.3f}) * (log10(4) - log10(2)) / ({F2_at_ha4:.3f} - {F2_at_ha2:.3f})")
    print(f"log(h/a)_design = {log_ha:.4f}")
    print(f"(h/a)_design = 10^{log_ha:.4f} = {h_a_design:.3f}\n")
    
    print("Finally, calculate the required pavement thickness h_design:")
    print(f"h_design = (h/a)_design * a_design = {h_a_design:.3f} * {a_design:.2f} mm")
    print(f"h_design = {h_design:.1f} mm")
    
    return h_design

if __name__ == '__main__':
    final_thickness = solve_pavement_thickness()
    print(f"\n<<<The required pavement thickness is {final_thickness:.1f} mm>>>")
    # For the final answer format requested
    # The submission system might not run the main block, so let's call the function and print the final format directly.
    # final_thickness = solve_pavement_thickness()
    # print(f"<<<{final_thickness:.1f}>>>")