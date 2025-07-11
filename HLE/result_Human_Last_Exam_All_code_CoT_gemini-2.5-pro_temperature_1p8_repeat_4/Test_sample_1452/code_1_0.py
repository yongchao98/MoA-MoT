def solve_equation():
    """
    This script identifies and prints the famous equation associated with Max Born.
    """
    
    # The scientist's name
    name = "Max Born"
    
    # The Born rule states that the probability P of finding a particle is the
    # square of the magnitude of its wave function (psi).
    # The number in this equation is the exponent.
    exponent = 2
    
    # Constructing the equation in LaTeX format.
    # We use f-string to insert the number into the equation.
    equation_latex = f"\\( P = |\\psi|^{{{exponent}}} \\)"
    
    # Print the final result
    print(f"Name: {name}")
    print(f"Equation: {equation_latex}")

solve_equation()