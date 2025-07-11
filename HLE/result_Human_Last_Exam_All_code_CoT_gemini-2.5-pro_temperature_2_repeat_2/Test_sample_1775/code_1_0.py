import numpy as np

def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the effective displacement field in a dual-gate FET.

    Args:
        Ctg (float): Top gate capacitance per unit area (in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Back gate capacitance per unit area (in F/m^2).
        Vbg (float): Back gate voltage (in V).
    """

    # The displacement field is the sum of the contributions from the top and bottom gates.
    # D = D_tg + D_bg
    # where D_tg = Ctg * Vtg and D_bg = Cbg * Vbg.
    # The final units are Coulombs per square meter (C/m^2).

    D_tg = Ctg * Vtg
    D_bg = Cbg * Vbg
    D_total = D_tg + D_bg

    print("The physical model treats the dual-gate transistor as two capacitors sharing a common grounded channel.")
    print("The total displacement field (D) is the sum of the contributions from the top gate (tg) and back gate (bg).")
    print("\nThe formula is: D = C_tg * V_tg + C_bg * V_bg\n")
    print("Substituting the given values into the equation:")
    # Using format specifiers to control the scientific notation and precision for clarity.
    print(f"D = ({Ctg:.2e} F/m^2) * ({Vtg:.2f} V) + ({Cbg:.2e} F/m^2) * ({Vbg:.2f} V)")
    print(f"D = {D_tg:.2e} C/m^2 + {D_bg:.2e} C/m^2")
    print(f"D = {D_total:.2e} C/m^2")
    
    # Return the value for potential further use, although the primary output is printing.
    return D_total

if __name__ == '__main__':
    # Example values for the parameters.
    # Capacitance values are typical for MOS devices.
    C_top_gate = 1.2e-6  # Capacitance per area in Farads per square meter (F/m^2)
    V_top_gate = 1.5     # Voltage in Volts (V)
    C_back_gate = 0.4e-6 # Capacitance per area in F/m^2
    V_back_gate = 5.0    # Voltage in V

    # Calculate and print the result
    final_D = calculate_displacement_field(C_top_gate, V_top_gate, C_back_gate, V_back_gate)
    # The final expression is embedded in the calculation function itself as requested.
    # For the final answer format, we derive it from the calculation.
    # D = 1.2e-6 * 1.5 + 0.4e-6 * 5.0 = 1.8e-6 + 2.0e-6 = 3.8e-6
    # Note: <<<...>>> should contain only the final value, not the equation.
    print(f"\n<<<{final_D:.1e}>>>")
