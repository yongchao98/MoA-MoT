def display_capacitance_formula():
    """
    This function prints the derived formula for the gate capacitance per unit area
    for the given quantum Hall system. The final formula expresses the capacitance
    in terms of the magnetic field (B), the base voltage (V_1), and fundamental
    constants e (elementary charge) and h (Planck's constant).
    """

    # The derived formula is C_g = 2 * e^2 * B / (h * V_1)
    
    print("Based on the physics of the Quantum Hall Effect, the gate capacitance per unit area (C_g) is calculated.")
    print("The final derived equation is:")
    
    # As requested, we explicitly output the numbers in the final equation.
    # The numerical coefficient in the numerator is 2.
    # The exponent for the elementary charge 'e' is 2.
    numerator_coefficient = 2
    charge_exponent = 2
    
    # Using an f-string to construct and print the formula.
    # The variables 'e', 'B', 'h', 'V_1' represent physical quantities.
    print(f"C_g = ({numerator_coefficient} * e^{charge_exponent} * B) / (h * V_1)")

# Execute the function to display the result
display_capacitance_formula()