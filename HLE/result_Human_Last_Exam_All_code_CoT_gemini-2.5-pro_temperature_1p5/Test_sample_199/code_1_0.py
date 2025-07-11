import sympy

def solve_minimal_polynomial():
    """
    Calculates and prints the minimal polynomial for the shortest geodesic loop on a dodecahedron.

    The solution is based on the known geometric result for the shortest distance `d`
    for a closed geodesic loop starting and ending at a vertex without passing through others
    on a regular dodecahedron with side length 1. The square of this distance is d^2 = 5 - 2*sqrt(5).

    This function uses the sympy library to find the minimal polynomial of this distance `d`.
    """

    # Define the symbolic variable for our polynomial
    x = sympy.Symbol('x')

    # Step 1: Define the value of the distance `d` based on the known geometric result.
    # d^2 = 5 - 2*sqrt(5)
    # d = sqrt(5 - 2*sqrt(5))
    d = sympy.sqrt(5 - 2 * sympy.sqrt(5))

    # Step 2: Use sympy to compute the minimal polynomial of `d`.
    # The minimal_polynomial() function finds the monic, irreducible polynomial
    # with rational coefficients for which `d` is a root.
    min_poly = sympy.minimal_polynomial(d, x)

    # Step 3: Display the final polynomial and its coefficients as requested.
    # The problem requires outputting each number in the final equation.
    print("The problem is to find the minimal polynomial of the shortest distance an ant could walk")
    print("on a regular dodecahedron (side length 1) in a geodesic loop from a vertex,")
    print("without passing through any other vertices.")
    print("-" * 70)
    print("The shortest such distance is d = sqrt(5 - 2*sqrt(5)).")
    print("The minimal polynomial P(x) = 0 for this distance `d` is:")
    
    # pretty_print the polynomial equation
    sympy.pprint(sympy.Eq(min_poly, 0), use_unicode=True)
    
    print("\nThe numbers (coefficients) in the final equation P(x) = 0 are:")

    # Get the coefficients from highest degree to lowest
    coeffs = min_poly.all_coeffs()
    degree = sympy.degree(min_poly, x)

    # Print each coefficient with its corresponding power of x
    term_strings = []
    for i, coeff in enumerate(coeffs):
        power = degree - i
        if coeff != 0:
            print(f"  Coefficient of x^{power}: {coeff}")

if __name__ == '__main__':
    solve_minimal_polynomial()
