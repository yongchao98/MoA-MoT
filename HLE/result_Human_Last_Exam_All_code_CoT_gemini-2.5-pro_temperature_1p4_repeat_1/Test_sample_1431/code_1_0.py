def solve_rope_challenge():
    """
    This function prints the derived formula for the force F required
    to lift the rope according to the problem's conditions.
    """
    
    # Symbolic representation of the variables
    m = 'm'  # mass of the rope
    l = 'l'  # length of the rope
    v = 'v'  # speed of the rope's center of mass
    g = 'g'  # acceleration due to gravity
    
    # Numerical coefficients from the derivation
    coeff_g_num = 2
    coeff_g_den = 3
    coeff_v = 2
    
    # Print the explanation and the final formula
    print("To conquer the mystical challenge, the monk must summon a precise constant force F.")
    print("Based on the principles of variable-mass dynamics, the formula for this force is:")
    print(f"F = ({coeff_g_num}/{coeff_g_den}) * {m} * {g} + {coeff_v} * {m} * {v}^2 / {l}")
    print("\nWhere:")
    print(f"  m = mass of the rope")
    print(f"  l = length of the rope")
    print(f"  v = speed of the rope's center of mass at the final moment")
    print(f"  g = acceleration due to gravity")

# Execute the function to display the answer
solve_rope_challenge()
