def calculate_displacement_field(Ctg, Vtg, Cbg, Vbg):
    """
    Calculates the displacement field in a dual-gate FET.

    Args:
        Ctg (float): Top gate capacitance per unit area (in F/m^2).
        Vtg (float): Top gate voltage (in V).
        Cbg (float): Back gate capacitance per unit area (in F/m^2).
        Vbg (float): Back gate voltage (in V).
    """

    # The displacement field D is the sum of the fields from the top and back gates.
    # D = D_top + D_bottom = Ctg * Vtg + Cbg * Vbg
    # The unit of displacement field is Coulombs per square meter (C/m^2).
    # Note: 1 Farad (F) = 1 Coulomb (C) / 1 Volt (V)

    D_total = Ctg * Vtg + Cbg * Vbg

    print("The formula for the displacement field (D) is:")
    print("D = C_tg * V_tg + C_bg * V_bg")
    print("\nUsing the provided values:")
    
    # Print the equation with the numbers substituted
    print(f"D = ({Ctg} F/m^2) * ({Vtg} V) + ({Cbg} F/m^2) * ({Vbg} V)")

    # Print the final result
    print(f"\nThe total displacement field through the transistor is:")
    print(f"D = {D_total} C/m^2")
    return D_total

if __name__ == '__main__':
    # --- Example Values ---
    # Typical gate capacitance is in the range of micro-Farads per cm^2.
    # We will use SI units (F/m^2). 1 uF/cm^2 = 1e-6 F / (1e-4 m^2) = 1e-2 F/m^2.

    # Let's assume C_tg = 2.0 uF/cm^2 and C_bg = 0.5 uF/cm^2.
    Ctg_si = 2.0E-2  # in F/m^2
    Vtg_si = 1.5     # in Volts
    Cbg_si = 0.5E-2  # in F/m^2
    Vbg_si = 5.0     # in Volts

    # Calculate and print the result
    final_answer = calculate_displacement_field(Ctg_si, Vtg_si, Cbg_si, Vbg_si)
    
    # The problem asks for the answer in a specific format at the end.
    # The calculated value based on the examples is 0.055.