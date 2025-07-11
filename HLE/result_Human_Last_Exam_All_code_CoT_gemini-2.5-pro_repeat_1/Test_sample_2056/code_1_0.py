import sympy

def solve_lk_n():
    """
    This function provides the symbolic derivation for l_k(n) and prints the result.
    """
    n, k = sympy.symbols('n k', real=True, positive=True)

    # The final derived formula for l_k(n)
    # l_k(n) = 1/2 * ln(n+1) - k^2 * (2 - 1/n)
    l_k_n_formula = sympy.Rational(1, 2) * sympy.log(n + 1) - k**2 * (2 - 1/n)
    
    # Expanding the formula
    l_k_n_expanded = sympy.expand(l_k_n_formula)

    print("The exact value of l_k(n) is:")
    sympy.pprint(l_k_n_expanded)

    # Outputting each numerical coefficient in the final equation as requested
    # The formula is of the form: c1*ln(c2*n + c3) + c4*k**2 + c5*k**2/n
    c1 = sympy.Rational(1, 2)
    c2 = 1
    c3 = 1
    c4 = -2
    c5 = 1
    
    print("\nCoefficients of the equation l_k(n) = c1*ln(c2*n + c3) + c4*k^2 + c5*k^2/n:")
    print(f"c1 (coefficient of the log term): {c1}")
    print(f"c2 (coefficient of n inside log): {c2}")
    print(f"c3 (constant inside log): {c3}")
    print(f"c4 (coefficient of k^2 term): {c4}")
    print(f"c5 (coefficient of k^2/n term): {c5}")

if __name__ == '__main__':
    solve_lk_n()