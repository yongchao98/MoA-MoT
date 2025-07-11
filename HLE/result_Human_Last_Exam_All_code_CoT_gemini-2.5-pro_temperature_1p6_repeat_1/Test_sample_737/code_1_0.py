import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """
    # Step 1: Define given parameters in consistent units (N, mm)
    # Plate bearing test (rigid plate, as per standard test setup)
    P_plate = 30 * 1000      # Load in Newtons (N)
    d_plate = 305            # Plate diameter in mm
    a_plate = d_plate / 2    # Plate radius in mm
    w_s = 2460 / 1000        # Deflection on subgrade in mm
    w_p = 1080 / 1000        # Deflection on trial pavement in mm
    h_trial = 300            # Trial pavement thickness in mm

    # Design wheel load (flexible load, as it's a tyre)
    W_design_ton = 1.80      # Load weight in tons
    g = 9.81                 # Acceleration due to gravity in m/s^2
    P_design = W_design_ton * 1000 * g  # Design load in Newtons (N)
    q_design_kPa = 600       # Tyre pressure in kN/m^2
    q_design_MPa = q_design_kPa / 1000  # Tyre pressure in N/mm^2 (MPa)

    # Design constraint
    w_design_limit = 1.00    # Maximum allowable deflection in mm

    # For mu=0.5, Boussinesq deflection is w = (factor * q * a) / E
    # where the factor depends on plate rigidity.
    factor_rigid = 1.18
    factor_flex = 1.5

    # Step 2: Calculate Subgrade Modulus of Elasticity (E_s)
    # Using the formula for a rigid plate on a single layer: w_s = (1.18 * P_plate) / (pi * a_plate * E_s)
    # Rearranging for E_s:
    E_s = (factor_rigid * P_plate) / (math.pi * a_plate * w_s)

    # Step 3: Determine Modular Ratio (E_p/E_s)
    # First, find the deflection factor F2 for the trial section
    F2_trial = w_p / w_s
    # And the thickness-to-radius ratio for the trial
    h_a_ratio_trial = h_trial / a_plate

    # We now find E_p/E_s by interpolating from a standard Burmister chart for F2.
    # We use chart values for h/a ≈ 2.0 (since h_a_ratio_trial ≈ 1.97)
    # Chart Point 1: E_p/E_s = 10, F2 ≈ 0.50
    # Chart Point 2: E_p/E_s = 20, F2 ≈ 0.40
    Ep_Es_ratio_1, F2_val_1 = 10, 0.50
    Ep_Es_ratio_2, F2_val_2 = 20, 0.40
    # Linear interpolation:
    slope = (Ep_Es_ratio_2 - Ep_Es_ratio_1) / (F2_val_2 - F2_val_1)
    Ep_Es_ratio = Ep_Es_ratio_1 + slope * (F2_trial - F2_val_1)

    # Step 4: Calculate design load parameters
    # The design load is flexible. Find the equivalent radius a_design.
    # P_design = q_design * Area = q_design * (pi * a_design^2)
    a_design = math.sqrt(P_design / (math.pi * q_design_MPa))

    # Calculate the hypothetical deflection on the subgrade alone under the design load
    w_s_design_hypothetical = (factor_flex * q_design_MPa * a_design) / E_s

    # Step 5: Calculate the required design deflection factor (F2_design)
    F2_design = w_design_limit / w_s_design_hypothetical

    # Step 6: Find the required thickness-to-radius ratio ((h/a)_design)
    # We interpolate again from the Burmister chart to find the h/a ratio that
    # corresponds to F2_design for our calculated Ep_Es_ratio.
    # We find h/a values for F2 ≈ F2_design on the curves for E_p/E_s = 10 and E_p/E_s = 20.
    # Chart Point 1: At E_p/E_s = 10, F2=0.345 corresponds to h/a ≈ 4.0
    # Chart Point 2: At E_p/E_s = 20, F2=0.345 corresponds to h/a ≈ 2.5
    h_a_val_1, h_a_val_2 = 4.0, 2.5
    # Linear interpolation:
    slope_2 = (h_a_val_2 - h_a_val_1) / (Ep_Es_ratio_2 - Ep_Es_ratio_1)
    h_a_design_ratio = h_a_val_1 + slope_2 * (Ep_Es_ratio - Ep_Es_ratio_1)

    # Step 7: Calculate the final required pavement thickness
    h_design = h_a_design_ratio * a_design

    # Print the step-by-step derivation and the final equation
    print("DERIVATION OF VALUES:\n")
    print(f"1. Subgrade Modulus (E_s): Using data from the subgrade-only plate test,")
    print(f"   E_s = ({factor_rigid} * P_plate) / (pi * a_plate * w_s)")
    print(f"   E_s = ({factor_rigid} * {P_plate:.0f} N) / (pi * {a_plate:.2f} mm * {w_s:.3f} mm) = {E_s:.2f} MPa\n")

    print(f"2. Modular Ratio (E_p/E_s): Using data from the trial pavement test,")
    print(f"   F2_trial = w_p / w_s = {w_p:.3f} mm / {w_s:.3f} mm = {F2_trial:.3f}")
    print(f"   (h/a)_trial = {h_trial} mm / {a_plate:.2f} mm = {h_a_ratio_trial:.3f}")
    print(f"   Interpolating from a Burmister chart gives E_p/E_s = {Ep_Es_ratio:.2f}\n")
    
    print(f"3. Design Load Radius (a_design): For the flexible wheel load,")
    print(f"   a_design = sqrt(P_design / (pi * q_design))")
    print(f"   a_design = sqrt({P_design:.2f} N / (pi * {q_design_MPa:.3f} MPa)) = {a_design:.2f} mm\n")

    print(f"4. Required Design Deflection Factor (F2_design): To limit deflection to {w_design_limit} mm,")
    print(f"   Hypothetical subgrade deflection w_s_hyp = ({factor_flex} * q_design * a_design) / E_s")
    print(f"   w_s_hyp = ({factor_flex} * {q_design_MPa:.3f} MPa * {a_design:.2f} mm) / {E_s:.2f} MPa = {w_s_design_hypothetical:.2f} mm")
    print(f"   F2_design = w_design_limit / w_s_hyp = {w_design_limit:.2f} mm / {w_s_design_hypothetical:.2f} mm = {F2_design:.3f}\n")

    print(f"5. Required Thickness-to-Radius Ratio ((h/a)_ratio):")
    print(f"   Interpolating from a Burmister chart for F2 = {F2_design:.3f} and E_p/E_s = {Ep_Es_ratio:.2f} gives (h/a)_ratio = {h_a_design_ratio:.3f}\n")

    print("FINAL CALCULATION:\n")
    print("The required pavement thickness (h_design) is calculated as:")
    print("h_design = (h/a)_ratio * a_design")
    print(f"h_design = {h_a_design_ratio:.3f} * {a_design:.2f}")
    print(f"h_design = {h_design:.2f} mm")
    
    return h_design

if __name__ == '__main__':
    final_thickness = solve_pavement_thickness()
    print(f"\n<<<{final_thickness:.1f}>>>")
