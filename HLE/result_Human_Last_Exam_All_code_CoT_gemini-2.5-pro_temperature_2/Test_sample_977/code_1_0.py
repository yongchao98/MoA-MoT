def print_potential_formulas():
    """
    Prints the derived symbolic formulas for the electric potential
    in the two regions between the grounded plates.
    """
    
    # Common denominator for the potential expressions
    denominator = "k * (epsilon_2 * cosh(k*a) * sinh(k*b) + epsilon_1 * sinh(k*a) * cosh(k*b))"
    
    # Numerator for the potential in Region 2 (0 < y < a)
    phi_2_numerator = "-sigma_0 * sinh(k*b) * sinh(k*(y - a)) * sin(k*x)"
    
    # Numerator for the potential in Region 1 (-b < y < 0)
    phi_1_numerator = "sigma_0 * sinh(k*a) * sinh(k*(y + b)) * sin(k*x)"
    
    # Construct the full expressions
    phi_2_expression = f"({phi_2_numerator}) / ({denominator})"
    phi_1_expression = f"({phi_1_numerator}) / ({denominator})"

    print("The derived electric potential Phi(x, y) is given by the following expressions:")
    print("-" * 70)
    
    print("\nIn the region 0 < y < a:")
    print(f"Phi(x, y) = {phi_2_expression}")
    
    print("\nIn the region -b < y < 0:")
    print(f"Phi(x, y) = {phi_1_expression}")
    
    print("-" * 70)
    print("\nThis result corresponds to answer choice A.")

if __name__ == '__main__':
    print_potential_formulas()
