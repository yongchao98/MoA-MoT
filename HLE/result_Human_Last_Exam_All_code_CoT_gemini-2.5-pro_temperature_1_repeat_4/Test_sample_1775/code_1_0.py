import sys

def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the displacement field in a dual-gate FET.

    Args:
        Ctg (float): Top gate capacitance per unit area (F/m^2).
        Vtg (float): Top gate voltage (V).
        Cbg (float): Back gate capacitance per unit area (F/m^2).
        Vbg (float): Back gate voltage (V).
    """
    # The displacement field D is equal to the total induced charge density
    # in the channel, which is the sum of charges from the top and back gates.
    # D = Ctg * Vtg + Cbg * Vbg

    charge_from_top_gate = Ctg * Vtg
    charge_from_back_gate = Cbg * Vbg
    displacement_field = charge_from_top_gate + charge_from_back_gate

    # Print the explanation and the final equation
    print("The displacement field (D) is calculated using the formula: D = Ctg * Vtg + Cbg * Vbg")
    print("\nUsing the provided values:")
    print(f"Ctg = {Ctg} F/m^2")
    print(f"Vtg = {Vtg} V")
    print(f"Cbg = {Cbg} F/m^2")
    print(f"Vbg = {Vbg} V")
    
    print("\nThe final equation with the numbers is:")
    print(f"D = ({Ctg} * {Vtg}) + ({Cbg} * {Vbg})")
    print(f"D = {charge_from_top_gate} + {charge_from_back_gate}")
    print(f"D = {displacement_field} C/m^2")
    
    # Output the final answer in the specified format
    # The 'file=sys.stderr' argument is used to separate the final answer
    # from the explanatory text, which is good practice for scripting.
    print(f"<<<{displacement_field}>>>", file=sys.stderr)


if __name__ == '__main__':
    # --- User-configurable values ---
    # Top gate capacitance per area in Farads per square meter (F/m^2)
    top_gate_capacitance_per_area = 1.5e-2

    # Top gate voltage in Volts (V)
    top_gate_voltage = 1.0

    # Back gate capacitance per area in Farads per square meter (F/m^2)
    back_gate_capacitance_per_area = 0.5e-2

    # Back gate voltage in Volts (V)
    back_gate_voltage = 2.0
    # ------------------------------------

    calculate_displacement_field(
        Ctg=top_gate_capacitance_per_area,
        Vtg=top_gate_voltage,
        Cbg=back_gate_capacitance_per_area,
        Vbg=back_gate_voltage
    )
