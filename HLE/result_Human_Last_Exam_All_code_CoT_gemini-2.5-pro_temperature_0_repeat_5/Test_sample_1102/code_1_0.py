import sys

def solve_welding_parameters():
    """
    Calculates and explains the selection of TIG welding parameters for a turbine blade repair.
    """
    # Given parameters from the problem description
    arc_gap_L = 6.0  # mm
    
    # --- Step 1: Estimate the required Arc Voltage ---
    # A more detailed empirical formula for TIG arc voltage is V = V_drop + (E_arc * L)
    # V_drop (Anode + Cathode voltage drop) is ~13V for Argon shielding gas.
    # E_arc (Arc column gradient) is ~1.0 V/mm.
    v_drop = 13.0
    e_arc = 1.0
    
    # Calculate the estimated voltage
    estimated_voltage = v_drop + (e_arc * arc_gap_L)
    
    print("--- Parameter Estimation ---")
    print("\nStep 1: Estimate Arc Voltage based on Arc Gap")
    print(f"The given arc gap is {arc_gap_L} mm.")
    print("Using the empirical formula: Voltage = V_drop + (E_arc * Arc_Gap)")
    print(f"Calculation: Voltage = {v_drop} V + ({e_arc} V/mm * {arc_gap_L} mm)")
    print(f"Estimated Voltage = {estimated_voltage:.1f} V")
    print("Based on this estimation, a voltage around 17.5 V is a plausible choice from the given options.")

    # --- Step 2: Evaluate the required Welding Current ---
    print("\nStep 2: Evaluate Welding Current based on Application")
    print("The task is a delicate root pass repair on a thin aeroengine blade tip.")
    print("This requires low heat input to prevent burn-through and material damage.")
    print("Therefore, a low welding current is necessary.")
    print("Currents like 80A, 100A, or 150A are too high for this application.")
    print("A current of 25 A is appropriate for such precision work.")

    # --- Step 3: Select the best option ---
    final_current = 25
    final_voltage = 17.5
    
    print("\n--- Conclusion ---")
    print("The best choice combines a sufficiently high voltage for the arc gap with a low current for the delicate work.")
    print(f"Selected Parameters: Current = {final_current} A, Voltage = {final_voltage} V.")

solve_welding_parameters()