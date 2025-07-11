import sys

def calculate_displacement_field(top_gate_capacitance, top_gate_voltage, back_gate_capacitance, back_gate_voltage):
    """
    Calculates the total displacement field in a dual-gate FET.

    Args:
        top_gate_capacitance (float): Capacitance per unit area of the top gate (in F/m^2).
        top_gate_voltage (float): Voltage applied to the top gate (in V).
        back_gate_capacitance (float): Capacitance per unit area of the back gate (in F/m^2).
        back_gate_voltage (float): Voltage applied to the back gate (in V).
    """
    # The displacement field is the sum of contributions from the top and back gates.
    # D = C_tg * V_tg + C_bg * V_bg
    # The unit of D will be (F/m^2) * V = (C/V/m^2) * V = C/m^2.
    
    total_displacement_field = top_gate_capacitance * top_gate_voltage + back_gate_capacitance * back_gate_voltage
    
    print("The formula for the total displacement field (D) is:")
    print("D = C_tg * V_tg + C_bg * V_bg\n")
    
    print("Plugging in the given values:")
    # Using 'g' format specifier to handle floating point and integer representation cleanly.
    print(f"D = {top_gate_capacitance:g} * {top_gate_voltage:g} + {back_gate_capacitance:g} * {back_gate_voltage:g}")
    
    term1 = top_gate_capacitance * top_gate_voltage
    term2 = back_gate_capacitance * back_gate_voltage
    
    print(f"D = {term1:g} + {term2:g}")
    print(f"D = {total_displacement_field:g} C/m^2")
    
    # Returning the final numerical answer as a string for the <<<>>> format
    return f"{total_displacement_field:g}"

if __name__ == '__main__':
    # --- User-defined variables ---
    # You can change these values to solve for your specific case.
    # C_tg: Top gate capacitance per unit area in Farads per square meter (F/m^2)
    C_tg = 1.2e-2
    # V_tg: Top gate voltage in Volts (V)
    V_tg = 1.0
    # C_bg: Back gate capacitance per unit area in Farads per square meter (F/m^2)
    C_bg = 0.5e-2
    # V_bg: Back gate voltage in Volts (V)
    V_bg = -2.0
    
    # Call the function to perform calculation and print the results
    final_answer = calculate_displacement_field(C_tg, V_tg, C_bg, V_bg)
    
    # The final answer is printed in the required format
    print(f"\n<<<{final_answer}>>>")
