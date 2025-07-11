def calculate_displacement_field(c_tg, v_tg, c_bg, v_bg):
    """
    Calculates the effective displacement field in a dual-gate FET.

    Args:
        c_tg (float): Top gate capacitance per unit area (F/m^2).
        v_tg (float): Top gate voltage (V).
        c_bg (float): Back gate capacitance per unit area (F/m^2).
        v_bg (float): Back gate voltage (V).

    Returns:
        float: The total effective displacement field (C/m^2).
    """

    # Displacement field from the top gate
    d_tg = c_tg * v_tg
    # Displacement field from the back gate
    d_bg = c_bg * v_bg
    # Total effective displacement field is the sum of the two contributions
    d_total = d_tg + d_bg

    print("The effective displacement field (D) through the transistor is the sum of the fields contributed by the top and back gates.")
    print("D = C_tg * V_tg + C_bg * V_bg")
    print(f"D = ({c_tg:.2e} F/m^2) * ({v_tg} V) + ({c_bg:.2e} F/m^2) * ({v_bg} V)")
    print(f"D = {d_tg:.2e} C/m^2 + {d_bg:.2e} C/m^2")
    print(f"Total Displacement Field (D) = {d_total:.2e} C/m^2")
    
    return d_total

if __name__ == '__main__':
    # --- User-defined variables ---
    # Top gate capacitance per area in Farads per square meter (F/m^2)
    C_tg = 1.7e-3
    # Back gate capacitance per area in Farads per square meter (F/m^2)
    C_bg = 4.5e-4
    # Top gate voltage in Volts (V)
    V_tg = 1.5
    # Back gate voltage in Volts (V)
    V_bg = 5.0
    
    # --- Calculation ---
    total_displacement_field = calculate_displacement_field(C_tg, V_tg, C_bg, V_bg)
    
    # Final answer formatting
    # print(f"<<<{total_displacement_field}>>>")
