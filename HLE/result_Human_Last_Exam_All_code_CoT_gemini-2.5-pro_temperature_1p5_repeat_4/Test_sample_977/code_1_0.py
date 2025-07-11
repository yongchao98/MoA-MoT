def generate_potential_formula():
    """
    This function generates and prints the mathematical formula for the electric potential
    in the region 0 < y < a, as derived from the physics principles.
    """
    
    # Define string representations of the variables for clarity in the final expression
    sigma_0 = "σ₀"
    k = "k"
    x = "x"
    y = "y"
    a = "a"
    b = "b"
    epsilon_1 = "ε₁"
    epsilon_2 = "ε₂"
    
    # Construct the numerator and denominator strings based on the derived formula
    numerator = f"-{sigma_0} sinh({k}{b}) sinh({k}({y} - {a})) sin({k}{x})"
    
    denominator_part1 = f"{epsilon_2} cosh({k}{a}) sinh({k}{b})"
    denominator_part2 = f"{epsilon_1} sinh({k}{a}) cosh({k}{b})"
    denominator = f"{k} [ {denominator_part1} + {denominator_part2} ]"
    
    # Print the full expression for the potential in the specified region
    print("The electric potential Φ(x, y) in the region 0 < y < a is:")
    print(f"Φ(x, y) = {numerator}")
    print(f"         {' ' * len('Φ(x, y) = ')}{'-' * len(numerator)}")
    print(f"         {' ' * len('Φ(x, y) = ')}{denominator}")

# Execute the function to display the formula
generate_potential_formula()