def display_potential_formula():
    """
    This function prints the derived mathematical expressions for the electric potential
    in the two regions.
    """
    
    # Common denominator for both expressions
    denominator = "k * [epsilon_2*cosh(k*a)*sinh(k*b) + epsilon_1*sinh(k*a)*cosh(k*b)]"
    
    # Numerator for the potential in region 2 (0 < y < a)
    numerator_2 = "-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    
    # Numerator for the potential in region 1 (-b < y < 0)
    numerator_1 = "sigma_0 * sinh(k*a) * sinh(k*(y + b)) * sin(k*x)"
    
    print("The electric potential Phi(x, y) is given by:")
    print("-" * 50)
    
    print("For the region 0 < y < a:")
    print(f"Phi_2(x, y) = ({numerator_2}) / ({denominator})")
    
    print("\nFor the region -b < y < 0:")
    print(f"Phi_1(x, y) = ({numerator_1}) / ({denominator})")
    
    print("-" * 50)
    print("\nThis result matches answer choice A.")

display_potential_formula()
