import sympy as sp

def solve_overlap_integral():
    """
    Presents the analytical solution for the 2s-2s overlap integral in H2+
    using the symbolic mathematics library Sympy.
    """
    # Define the symbols for internuclear distance (R) and effective nuclear charge (zeta)
    # They are defined as positive real numbers.
    R, zeta = sp.symbols('R zeta', positive=True, real=True)

    print("This script presents the final analytical expression for the overlap integral (S) of two 2s orbitals.")
    print("The expression is given in terms of the internuclear distance R and the effective nuclear charge zeta.")
    print("-" * 70)

    # Define the dimensionless variable rho for simplicity
    rho = zeta * R / 2

    # The final equation is composed of an exponential part and a polynomial part.
    # Let's define each term of the polynomial.
    exponential_term = sp.exp(-rho)
    term_1 = 1
    term_2 = rho
    term_3 = rho**2 / 3
    term_4 = rho**4 / 15

    # Assemble the full equation
    S_final = exponential_term * (term_1 + term_2 + term_3 + term_4)

    # Print the equation in a structured way, showing each component.
    print("The final equation for the overlap integral S is:")
    print("\nS = (Exponential Term) * (Term 1 + Term 2 + Term 3 + Term 4)\n")
    print("Where each component is defined as follows:\n")

    print("Exponential Term:")
    sp.pprint(exponential_term, use_unicode=True)
    print("\nTerm 1:")
    sp.pprint(term_1, use_unicode=True)
    print("\nTerm 2:")
    sp.pprint(term_2, use_unicode=True)
    print("\nTerm 3:")
    sp.pprint(term_3, use_unicode=True)
    print("\nTerm 4:")
    sp.pprint(term_4, use_unicode=True)

    print("-" * 70)
    print("The full expression for S(R, zeta) is:")
    sp.pprint(S_final, use_unicode=True)


if __name__ == "__main__":
    solve_overlap_integral()
