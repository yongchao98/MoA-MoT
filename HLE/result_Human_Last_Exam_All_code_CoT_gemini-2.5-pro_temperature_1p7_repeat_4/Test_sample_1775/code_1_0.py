import sys

def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the net displacement field in a dual-gate transistor.

    Args:
        Ctg (float): Top gate capacitance per area (in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Back gate capacitance per area (in F/m^2).
        Vbg (float): Back gate voltage (in V).
    """

    # The total displacement field is the superposition of the fields from the
    # back and top gates. A positive field is defined as pointing from the
    # back gate towards the top gate.
    # D = D_from_back_gate + D_from_top_gate
    # D = (C_bg * V_bg) - (C_tg * V_tg)
    
    # Calculate components
    D_bg_term = Cbg * Vbg
    D_tg_term = Ctg * Vtg
    
    # Calculate total displacement field
    D_total = D_bg_term - D_tg_term

    # Output the explanation and step-by-step calculation
    print("The net displacement field (D) is calculated using the formula: D = C_bg * V_bg - C_tg * V_tg")
    print("\n--- Calculation Steps ---")
    print(f"Given values:")
    print(f"  Top Gate Capacitance per Area (C_tg): {Ctg:.2e} F/m^2")
    print(f"  Top Gate Voltage (V_tg): {Vtg:.2f} V")
    print(f"  Back Gate Capacitance per Area (C_bg): {Cbg:.2e} F/m^2")
    print(f"  Back Gate Voltage (V_bg): {Vbg:.2f} V\n")

    print("Equation with values:")
    print(f"D = ({Cbg:.2e} F/m^2) * ({Vbg:.2f} V) - ({Ctg:.2e} F/m^2) * ({Vtg:.2f} V)")
    
    print("\nResulting in:")
    print(f"D = {D_bg_term:.2e} C/m^2 - {D_tg_term:.2e} C/m^2")

    print("\n--- Final Answer ---")
    print(f"The net displacement field through the transistor is: {D_total:.2e} C/m^2")
    
if __name__ == '__main__':
    # --- Example Values ---
    # Top gate capacitance per area (e.g., 1.5 uF/cm^2 = 1.5e-2 F/m^2)
    Ctg_val = 1.5e-2
    # Top gate voltage
    Vtg_val = 1.0
    # Back gate capacitance per area (e.g., 0.5 uF/cm^2 = 0.5e-2 F/m^2)
    Cbg_val = 0.5e-2
    # Back gate voltage
    Vbg_val = -2.0

    calculate_displacement_field(Ctg=Ctg_val, Vtg=Vtg_val, Cbg=Cbg_val, Vbg=Vbg_val)