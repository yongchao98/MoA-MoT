import math

def solve_pavement_thickness():
    """
    Calculates the required pavement thickness using Burmister's two-layer theory.
    """

    # --- 1. Input Parameters ---
    # Plate Bearing Test Data
    d_plate = 305.0  # mm
    P_plate = 30000.0  # N (30 kN)
    delta_s = 2460.0 / 1000  # mm (deflection on subgrade)
    delta_p = 1080.0 / 1000  # mm (deflection on trial pavement)
    h_trial = 300.0  # mm (trial pavement thickness)

    # Design Wheel Load Data
    W_design_ton = 1.80  # ton
    g = 9.81 # m/s^2
    P_design = W_design_ton * 1000 * g # N
    q_design_kpa = 600.0  # kN/m^2
    
    # Design Constraint
    delta_design_limit = 1.00  # mm

    # --- Calculations ---
    
    # Convert design pressure to N/mm^2 (MPa)
    q_design = q_design_kpa / 1000 # N/mm^2

    # Radius of the loading plate
    a_plate = d_plate / 2.0
    
    # Pressure under the loading plate
    q_plate = P_plate / (math.pi * a_plate**2)

    print("--- Step 1: Analyze Plate Bearing Test ---")
    print(f"Plate Radius (a_plate): {a_plate:.2f} mm")
    print(f"Plate Pressure (q_plate): {q_plate:.4f} N/mm^2")
    
    # --- 2. Calculate Subgrade Modulus (E_s) ---
    # Using flexible plate formula: delta_s = 1.5 * q_plate * a_plate / E_s
    E_s = (1.5 * q_plate * a_plate) / delta_s
    print(f"\n--- Step 2: Calculate Subgrade Modulus (E_s) ---")
    print(f"Subgrade Modulus (E_s): {E_s:.2f} N/mm^2 (MPa)")

    # --- 3. Determine Modular Ratio (E_p/E_s) ---
    # F2 for the trial section is the ratio of deflections
    F2_trial = delta_p / delta_s
    h_a_trial_ratio = h_trial / a_plate
    
    # From Burmister's charts, for F2 ≈ 0.439 and h/a ≈ 1.97, the modular ratio E_p/E_s is ~10.
    Ep_Es_ratio = 10.0
    
    print(f"\n--- Step 3: Determine Modular Ratio (E_p/E_s) ---")
    print(f"Trial Deflection Factor (F2_trial): {F2_trial:.4f}")
    print(f"Trial Geometry Ratio (h/a_trial): {h_a_trial_ratio:.4f}")
    print(f"Assumed Modular Ratio (E_p/E_s) from chart: {Ep_Es_ratio}")

    # --- 4. Calculate Design Load Parameters ---
    # Area = P_design / q_design; Area = pi * a_design^2
    a_design = math.sqrt(P_design / (math.pi * q_design))
    print(f"\n--- Step 4: Calculate Design Load Parameters ---")
    print(f"Design Load (P_design): {P_design:.2f} N")
    print(f"Design Load Radius (a_design): {a_design:.2f} mm")

    # --- 5. Calculate Required Deflection Factor (F2_design) ---
    # delta_design = (1.5 * q_design * a_design / E_s) * F2_design
    # F2_design = delta_design * E_s / (1.5 * q_design * a_design)
    F2_design = (delta_design_limit * E_s) / (1.5 * q_design * a_design)
    print(f"\n--- Step 5: Calculate Required Deflection Factor (F2_design) ---")
    print(f"Required Deflection Factor (F2_design): {F2_design:.4f}")

    # --- 6. Determine Required Pavement Thickness (h_design) ---
    def get_h_a_from_F2(F2_target, Ep_Es_ratio_in):
        """
        Simulates lookup from Burmister's F2 chart using linear interpolation.
        Finds h/a for a given F2 on the E_p/E_s = 10 curve.
        """
        chart_data = {
            10: [ (2.0, 0.43), (1.5, 0.53) ] # (h/a, F2) points for E_p/E_s = 10
        }
        points = chart_data.get(Ep_Es_ratio_in)
        if not points: return None
        
        p1, p2 = points[0], points[1] # p1=(h/a_1, F2_1), p2=(h/a_2, F2_2)
        h_a1, F2_1 = p1
        h_a2, F2_2 = p2
        
        # Linear interpolation: h_a = h_a1 + (h_a2 - h_a1) * (F2_target - F2_1) / (F2_2 - F2_1)
        h_a = h_a1 + (h_a2 - h_a1) * (F2_target - F2_1) / (F2_2 - F2_1)
        return h_a

    h_a_design_ratio = get_h_a_from_F2(F2_design, Ep_Es_ratio)
    h_design = h_a_design_ratio * a_design
    
    print(f"\n--- Step 6: Determine Required Pavement Thickness (h_design) ---")
    print(f"Required Geometry Ratio (h/a_design) from chart simulation: {h_a_design_ratio:.4f}")
    
    print("\n--- Final Calculation ---")
    print(f"The required pavement thickness is calculated as: h_design = (h/a)_design * a_design")
    print(f"h_design = {h_a_design_ratio:.4f} * {a_design:.2f}")
    print(f"Required Pavement Thickness (h_design) = {h_design:.2f} mm")
    
    return h_design

# Run the solver and print the final answer in the required format
final_thickness = solve_pavement_thickness()
print(f"\n<<<The required pavement thickness is {final_thickness:.2f} mm>>>")
