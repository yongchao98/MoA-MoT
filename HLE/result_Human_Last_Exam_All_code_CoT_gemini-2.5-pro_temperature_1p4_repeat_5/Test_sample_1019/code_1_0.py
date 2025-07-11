import sympy

def solve_hamiltonicity_threshold():
    """
    This function determines and displays the d-threshold for Hamiltonicity
    for a graph H_n U G(n, p) based on the given conditions.
    """
    
    # Define symbolic variables for n and eta
    n, eta = sympy.symbols('n eta')

    # The problem states the minimum degree d of the graph H_n is given by the equation:
    d = n/2 - eta
    
    print("The minimum degree d is defined by the equation:")
    sympy.pprint(sympy.Eq(sympy.Symbol('d'), d))
    print("-" * 20)
    
    # The d-threshold for Hamiltonicity is based on a result by Ben-Shimon, Krivelevich, and Sudakov.
    # Their (corrected) theorem states that for d = (1/2 - alpha)*n, the threshold p is
    # Theta(ln(1/alpha) / (alpha * n)).
    
    # From d = n/2 - eta, we can identify alpha = eta/n.
    # Substituting alpha into the formula gives:
    # p = Theta(ln(1/(eta/n)) / ((eta/n)*n))
    #   = Theta(ln(n/eta) / eta)
    
    # We will represent the threshold p by the core of this expression.
    # In this context, 'log' denotes the natural logarithm.
    p_threshold = sympy.log(n / eta) / eta

    print("The d-threshold for Hamiltonicity, p, is given by the expression:")
    sympy.pprint(sympy.Eq(sympy.Symbol('p'), p_threshold))
    
    # The problem specifies the range for eta. We can state it for context.
    print("\nThis threshold is valid for eta in the range [1/2, n/64].")


if __name__ == '__main__':
    solve_hamiltonicity_threshold()
