def calculate_displacement_field():
    """
    Calculates the displacement field in a dual-gate FET with a grounded channel.
    """
    # --- Input Variables (with example values) ---
    # Top gate voltage in Volts (V)
    V_tg = 1.5
    # Back gate voltage in Volts (V)
    V_bg = 0.8
    # Top gate capacitance per unit area in Farads per square meter (F/m^2)
    C_tg = 2.5e-3
    # Back gate capacitance per unit area in Farads per square meter (F/m^2)
    C_bg = 1.2e-3

    # --- Explanation ---
    print("In a dual-gate FET with a grounded channel, the total displacement field (D) is the sum of the contributions from the top and back gates.")
    print("The displacement field from each gate is given by D_gate = C_gate * V_gate, where C is the capacitance per area and V is the gate voltage.")
    print("The final formula is: D = C_tg * V_tg + C_bg * V_bg")
    print("-" * 40)

    # --- Calculation ---
    # The calculation is based on the principle of linear superposition.
    D = (C_tg * V_tg) + (C_bg * V_bg)

    # --- Output Results ---
    print("Using the provided example values:")
    print(f"C_tg = {C_tg} F/m^2")
    print(f"V_tg = {V_tg} V")
    print(f"C_bg = {C_bg} F/m^2")
    print(f"V_bg = {V_bg} V")
    print("\nThe final equation with these values is:")
    # Here we output each number in the final equation as requested
    print(f"D = ({C_tg}) * ({V_tg}) + ({C_bg}) * ({V_bg})")

    # The result is in Coulombs per square meter (C/m^2)
    print(f"\nThe resulting displacement field D is: {D} C/m^2")

# Execute the function
calculate_displacement_field()