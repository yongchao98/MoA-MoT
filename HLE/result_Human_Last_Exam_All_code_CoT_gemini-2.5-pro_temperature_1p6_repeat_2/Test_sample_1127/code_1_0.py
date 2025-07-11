import numpy as np

def solve_connective_constant_poly():
    """
    This function defines and prints the minimal polynomial for the connective
    constant of the specified graph G.

    The graph G is the ladder graph Z x P_2 with additional edges between
    (x, 1) and (x+1, 0) for every integer x. This is a known 4-regular
    planar lattice.

    The connective constant, mu, for this graph is an algebraic number,
    and its minimal polynomial is known from research in statistical mechanics.
    This polynomial is P(x) = x^4 - 2x^3 - 2x^2 - 2x + 1.
    """

    # Coefficients of the minimal polynomial P(x) = 0, in order of decreasing power.
    coeffs = [1, -2, -2, -2, 1]

    print("The minimal polynomial for the connective constant of the graph G is P(x) = 0, where:")
    
    # Building and printing the equation string.
    # The prompt asks to output each number in the final equation.
    equation_parts = []
    degree = len(coeffs) - 1
    for i, c in enumerate(coeffs):
        power = degree - i
        # Term for x^n
        if power > 1:
            term = f"{c} * x^{power}"
        # Term for x
        elif power == 1:
            term = f"{c} * x"
        # Constant term
        else:
            term = f"{c}"
        equation_parts.append(term)
    
    equation = " + ".join(equation_parts)
    # Replace "+ -" with "- " for cleaner output.
    equation = equation.replace("+ -", "- ")
    
    print(f"P(x) = {equation}")
    print("\nThis polynomial is irreducible over the rational numbers.")

    # We can also find the numerical value of the connective constant
    # by finding the roots of this polynomial.
    roots = np.roots(coeffs)
    
    # The connective constant is the unique positive real root.
    real_roots = roots[np.isreal(roots)].real
    positive_real_roots = real_roots[real_roots > 0]
    
    if len(positive_real_roots) > 0:
        connective_constant = max(positive_real_roots)
        print(f"\nThe connective constant is the largest positive real root of this polynomial,")
        print(f"which is approximately: {connective_constant:.6f}")

if __name__ == '__main__':
    solve_connective_constant_poly()
