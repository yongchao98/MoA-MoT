def print_scientist_equation():
    """
    Prints the name and famous equation for Max Born.
    """
    # Define the scientist's name. The prompt used quotes.
    scientist_name = '"Max Born"'

    # Define the equation associated with the scientist.
    # The Born rule states that the probability density is equal to the square of the absolute value of the wave function (Psi).
    # The individual components of the equation symbol are the vertical bar '|', the symbol '\Psi', and the exponent '^2'.
    equation_latex = "|\\Psi|^2"

    # Print the name and the equation in the specified format.
    print(f"Name: {scientist_name}")
    print(f"Equation: \\( {equation_latex} \\)")

print_scientist_equation()