def solve_displacement_field():
    """
    This function derives and prints the equation for the displacement field
    in a dual-gate FET.
    """
    # Define symbolic variables for the equation components.
    # C_tg: Top gate capacitance per unit area
    # V_tg: Top gate voltage
    # C_bg: Back gate capacitance per unit area
    # V_bg: Back gate voltage
    C_tg = "C_tg"
    V_tg = "V_tg"
    C_bg = "C_bg"
    V_bg = "V_bg"

    # The displacement field (D) is equal to the total surface charge density
    # induced in the transistor's channel.
    # The total charge is the sum of the charge induced by the top gate and the back gate.
    # Charge from top gate = C_tg * V_tg
    # Charge from back gate = C_bg * V_bg
    
    # We construct the final equation as a string.
    # The equation is D = C_tg*V_tg + C_bg*V_bg
    
    print("The equation for the displacement field (D) is:")
    print(f"D = {C_tg} * {V_tg} + {C_bg} * {V_bg}")

# Execute the function to get the result.
solve_displacement_field()
