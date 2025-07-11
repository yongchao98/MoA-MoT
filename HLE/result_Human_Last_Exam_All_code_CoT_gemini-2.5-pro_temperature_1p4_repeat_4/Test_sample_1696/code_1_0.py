def print_final_formula():
    """
    Prints the derived closed-form formula for P(n) and its constituent numbers.
    """
    
    # The derived closed formula is P(n) = (2 * pi)^((n * (n - 1)) / 4) / sqrt(n!)
    
    print("The derived closed formula for P(n) is expressed as:")
    print("P(n) = (2 * pi) ^ ( (n * (n - 1)) / 4 ) / ( (n!) ^ (1 / 2) )")
    
    print("\nAs requested, here are the constant numbers that form this equation:")
    
    # Numbers from the numerator part: (2 * pi)^(...)
    num_base = 2
    print(f"In the base of the numerator (2 * pi), the number is: {num_base}")
    
    # Numbers from the exponent of the numerator: (n * (n-1)) / 4
    exp_num_coeff = 1
    exp_den = 4
    print(f"In the exponent of the numerator ( (n * (n - {exp_num_coeff})) / {exp_den} ), the numbers are: {exp_num_coeff}, {exp_den}")
    
    # Numbers from the exponent of the denominator: (n!)^(1/2)
    den_exp_num = 1
    den_exp_den = 2
    print(f"In the exponent of the denominator (1 / 2), the numbers are: {den_exp_num}, {den_exp_den}")

print_final_formula()
