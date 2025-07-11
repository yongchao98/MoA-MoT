def print_magnetization_curve():
    """
    This function prints the derived analytical expression for the initial
    magnetization curve M(H) of a superconducting slab in a parallel field.
    """
    
    # Define the symbols used in the equations
    # M: Magnetization
    # H: Applied magnetic field
    # Jc: Critical current density
    # a: Half-thickness of the slab
    
    # The full penetration field, H_p, is the field at which the shielding
    # currents reach the center of the slab. H_p = Jc * a.
    
    print("The analytical expression for the initial magnetization curve M(H) is derived from the critical-state model.")
    print("The expression is piecewise, depending on whether the applied field H is below or above the full penetration field (Jc * a).\n")
    
    # --- Regime 1: Partial Field Penetration ---
    print("1. For 0 <= H <= Jc * a (Partial Penetration):")
    # In this regime, the magnetization is a quadratic function of H.
    # The equation shows an initial diamagnetic response (M = -H) corrected by a positive term.
    equation1 = "M(H) = -H + H^2 / (2 * a * Jc)"
    print(f"   {equation1}\n")
    
    # --- Regime 2: Full Field Penetration ---
    print("2. For H > Jc * a (Full Penetration):")
    # In this regime, the superconductor is fully saturated with shielding currents,
    # and the magnetization reaches a constant, maximum negative value.
    equation2 = "M(H) = -(Jc * a) / 2"
    print(f"   {equation2}\n")

# Execute the function to print the result
print_magnetization_curve()