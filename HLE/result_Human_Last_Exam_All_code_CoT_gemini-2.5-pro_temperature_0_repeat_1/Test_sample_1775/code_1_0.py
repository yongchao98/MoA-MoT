def calculate_displacement_field():
    """
    Calculates the net displacement field in a dual-gate FET.

    The total displacement field (D) terminating on the grounded channel is the
    linear superposition of the fields from the top and back gates.
    The formula is D = C_tg * V_tg + C_bg * V_bg.
    """

    # --- User-defined variables ---
    # C_tg: Top gate capacitance per unit area (in Farads per square meter, F/m^2)
    Ctg = 2.5e-3  # Example value: 2.5 mF/m^2

    # V_tg: Top gate voltage (in Volts, V)
    Vtg = 1.5  # Example value: 1.5 V

    # C_bg: Back gate capacitance per unit area (in Farads per square meter, F/m^2)
    Cbg = 0.8e-3  # Example value: 0.8 mF/m^2

    # V_bg: Back gate voltage (in Volts, V)
    Vbg = -5.0  # Example value: -5.0 V

    # --- Calculation ---
    # The total displacement field is the sum of the fields from each gate.
    # This quantity is equivalent to the total charge density induced in the channel.
    D_total = Ctg * Vtg + Cbg * Vbg

    # --- Output ---
    # Print the final equation with all the numbers
    print("The displacement field (D) is calculated as the sum of the contributions from the top and back gates.")
    print("Formula: D = C_tg * V_tg + C_bg * V_bg")
    print(f"D = {Ctg} F/m^2 * {Vtg} V + {Cbg} F/m^2 * {Vbg} V")
    print(f"D = {Ctg * Vtg} C/m^2 + {Cbg * Vbg} C/m^2")
    print(f"D = {D_total} C/m^2")

    # The final answer in the required format
    print(f"\n<<<{D_total}>>>")

calculate_displacement_field()