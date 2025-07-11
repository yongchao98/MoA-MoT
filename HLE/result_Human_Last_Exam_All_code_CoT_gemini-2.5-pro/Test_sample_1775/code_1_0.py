def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the net displacement field in a dual-gate FET.

    Args:
        Ctg (float): Top gate capacitance per unit area (in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Back gate capacitance per unit area (in F/m^2).
        Vbg (float): Back gate voltage (in V).
    """

    # The displacement field is the superposition of the fields from the top and back gates.
    # D = D_bg - D_tg
    # Let's define the positive direction as upward (from back gate to top gate).
    # The field from the top gate points down, so it's negative in our system.
    # The field from the back gate points up, so it's positive.
    # D_tg = Ctg * Vtg
    # D_bg = Cbg * Vbg
    # Net D = D_bg - D_tg
    
    D_net = Cbg * Vbg - Ctg * Vtg

    print("--- Calculating Displacement Field (D) ---")
    print(f"Top Gate Capacitance per Area (Ctg): {Ctg} F/m^2")
    print(f"Top Gate Voltage (Vtg): {Vtg} V")
    print(f"Back Gate Capacitance per Area (Cbg): {Cbg} F/m^2")
    print(f"Back Gate Voltage (Vbg): {Vbg} V")
    print("\nFormula: D = (Cbg * Vbg) - (Ctg * Vtg)")
    
    # As requested, outputting each number in the final equation
    print(f"Equation with values: D = ({Cbg} * {Vbg}) - ({Ctg} * {Vtg})")
    
    print(f"\nResult:")
    print(f"The net displacement field (D) is: {D_net:.4g} C/m^2")
    if D_net > 0:
        print("(Positive sign indicates the field is directed from the back gate to the top gate)")
    elif D_net < 0:
        print("(Negative sign indicates the field is directed from the top gate to the back gate)")
    else:
        print("(The net displacement field is zero)")

# --- Example Values ---
# Note: 1 uF/cm^2 = 1e-6 F / (1e-2 m)^2 = 1e-6 F / 1e-4 m^2 = 0.01 F/m^2
C_tg_example = 1.5e-2  # Corresponds to 1.5 uF/cm^2
V_tg_example = 2.0     # 2.0 Volts
C_bg_example = 0.5e-2  # Corresponds to 0.5 uF/cm^2
V_bg_example = -10.0   # -10.0 Volts

# Execute the calculation with example values
calculate_displacement_field(C_tg_example, V_tg_example, C_bg_example, V_bg_example)

# The final answer is the formula for the displacement field.
# Since the problem is symbolic, the formula is the main result.
# Let's represent the final answer in symbolic form.
final_answer = "D = C_bg*V_bg - C_tg*V_tg"
# This format is just for the thought process; the code above provides the full answer.