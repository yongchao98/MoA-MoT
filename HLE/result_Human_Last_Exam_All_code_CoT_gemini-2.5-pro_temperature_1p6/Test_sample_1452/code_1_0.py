def print_max_born_equation():
    """
    This function prints the name "Max Born" and his famous equation,
    the normalization condition derived from the Born rule.
    """
    name = "Max Born"
    
    # The equation is the normalization condition for the wave function.
    # It states that the integral of the probability density (|Ψ|^2) over all space (dτ) is 1.
    # The numbers in the equation are the exponent '2' and the result '1'.
    equation_string = r"\int |\Psi|^2 d\tau = 1"
    
    print(f"Name: {name}")
    print(f"Equation: \\({equation_string}\\)")

print_max_born_equation()