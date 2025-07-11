def print_force_equation():
    """
    This function prints the derived symbolic equation for the force per unit area on the conductor.
    """
    
    # Define the components of the equation as strings
    term_prefactor = "(\u03BC\u2080 / 2)"  # mu_0 / 2
    term_numerator = "K\u2080\u00B2 sin\u00B2(ay)"  # K_0^2 sin^2(ay)
    term_denominator = "[\n\t\t cosh(ad) + (\u03BC\u2080/\u03BC) sinh(ad) \n\t    ]^2" # [cosh(ad) + (mu_0/mu)sinh(ad)]^2
    term_direction = "\u00EE\u2093" # i_x unit vector
    
    # Assemble the final equation string
    # Using unicode for clarity of mathematical symbols
    
    line1 = f"  ->        {term_numerator}"
    line2 =  "  f / area = " + term_prefactor + " *  " + "-" * (len(term_numerator) + 2) + " * " + term_direction
    line3 = f"            {term_denominator}"
    
    print("The derived expression for the force per unit area is:")
    print("")
    print(line1)
    print(line2)
    print(line3)
    print("")

print_force_equation()