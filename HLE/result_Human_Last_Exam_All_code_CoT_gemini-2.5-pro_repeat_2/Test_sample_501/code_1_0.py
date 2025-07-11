import sympy

def polymer_force_law():
    """
    This function uses the sympy library to derive and display the force law
    for a thermally isolated freely jointed polymer chain.
    """
    # Define the symbols for the variables in the problem
    E0 = sympy.Symbol('E(0)', positive=True, real=True)
    x = sympy.Symbol('x', positive=True, real=True)
    n = sympy.Symbol('n', integer=True, positive=True)
    l = sympy.Symbol('l', positive=True, real=True)
    
    # The problem asks for the force law, which is the tension 'tau'.
    # Based on the derivation, the force law is:
    # tau = (2 * E(0) * x / (n**2 * l**2)) * exp(x**2 / (n**2 * l**2))
    
    # We construct this expression using sympy symbols
    numerator_factor = 2 * E0 * x
    denominator = n**2 * l**2
    
    exponential_term = sympy.exp(x**2 / denominator)
    
    force_tau = (numerator_factor / denominator) * exponential_term
    
    # Print the final force law
    print("The force law tau(x) is given by:")
    sympy.pprint(force_tau)
    
    # The problem asks to output each number in the final equation.
    # The numbers are 2 (in the numerator) and 2 (in the exponents).
    # Sympy's pretty print already does this.
    # We can also print the equation as a string for clarity.
    print("\nIn string format:")
    print(f"τ(x) = {force_tau}")

if __name__ == '__main__':
    polymer_force_law()
