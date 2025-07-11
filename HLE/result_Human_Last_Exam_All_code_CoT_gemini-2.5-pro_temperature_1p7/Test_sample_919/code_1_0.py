def print_final_equation():
    """
    Prints the derived formula for the force per unit area on the conductor.
    The formula is broken down into its constituent parts for clarity.
    """
    
    # The different symbolic parts of the final equation
    coefficient = "- (mu_0 / 2)"
    numerator_term = "K_0**2 * sin(a*y)**2"
    denominator_term = "[cosh(a*d) + (mu_0/mu) * sinh(a*d)]**2"
    direction_vector = "i_x"
    
    print("The final equation for the force per unit area (f/area) is constructed from the following parts:")
    print("-" * 70)
    
    # Printing each component as requested
    print(f"1. The constant coefficient:")
    print(f"   {coefficient}")
    
    print(f"\n2. The numerator term, dependent on the current sheet amplitude and position:")
    print(f"   {numerator_term}")

    print(f"\n3. The denominator term, representing the shielding effect of the material and air gap:")
    print(f"   {denominator_term}")

    print(f"\n4. The direction of the force vector:")
    print(f"   {direction_vector}")
    
    print("-" * 70)
    print("\nCombining these parts, the final equation is:")
    print(f"f/area = {coefficient} * ( {numerator_term} ) / ( {denominator_term} ) * {direction_vector}")

if __name__ == '__main__':
    print_final_equation()