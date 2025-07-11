import sympy

def find_hamiltonicity_threshold():
    """
    This function calculates and prints the formula for the d-threshold for Hamiltonicity
    based on the given parameters.
    """
    # Define the symbols for n and eta
    n = sympy.Symbol('n')
    eta = sympy.Symbol('eta')

    # The minimum degree d is given as a function of n and eta
    # We use sympy.Rational for the fraction to maintain precision
    d = sympy.Rational(1, 2) * n - eta

    # According to theoretical results (Krivelevich, 2005), the threshold p is
    # proportional to 1 divided by (n - d - 1) * n.
    # We construct this denominator.
    denominator = (n - d - 1) * n

    # To fulfill the requirement of showing each number in the final equation,
    # we expand the expression for the denominator.
    expanded_denominator = sympy.expand(denominator)
    
    # We now print the results in a clear format.
    # The threshold p is proportional to 1 over this denominator.
    print("The d-threshold for Hamiltonicity, denoted by p, is proportional to 1/D, where D is given by the formula:")

    # Printing the expanded formula highlights all the numerical coefficients.
    # The coefficients are 1/2 for n**2, 1 for eta*n, and -1 for n.
    print(f"D = ({sympy.S(1)/2})*n**2 + (1)*n*eta + (-1)*n")

    # For a more compact representation, we can also show the factored form.
    simplified_denominator = sympy.factor(expanded_denominator)
    print("\nIn a more compact, factored form, the denominator is:")
    print(f"D = {simplified_denominator}")


if __name__ == '__main__':
    find_hamiltonicity_threshold()
