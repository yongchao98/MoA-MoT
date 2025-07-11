import sympy

def solve_seifert_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles of the 9_23 knot
    using its HOMFLY polynomial.
    """
    # Step 1: Define the HOMFLY polynomial for the 9_23 knot.
    # This is a known invariant, which we define as a symbolic expression.
    a, z = sympy.symbols('a z')
    # The polynomial is taken from standard knot theory tables (e.g., Knot Atlas).
    homfly_poly_str = "-a**6/z - a**6*z + 2*a**4/z - a**4 + a**4*z - a**2/z + 2*a**2 - a**2*z"
    homfly_poly = sympy.sympify(homfly_poly_str)

    print(f"To find a lower bound for the minimum number of Seifert circles of the 9_23 knot, we use its HOMFLY polynomial.")
    print(f"The HOMFLY polynomial for 9_23 is P(a, z) = {homfly_poly_str}")
    print("\nWe use the inequality: s(K) >= span_a(P) / 2 + 1")
    print("where s(K) is the minimum number of Seifert circles and span_a(P) is the span of the polynomial with respect to variable 'a'.")

    # Step 2: Calculate the span of the polynomial in the variable 'a'.
    # We can treat the expression as a polynomial in 'a' to find the degrees.
    poly_in_a = sympy.Poly(homfly_poly, a)

    # Find the maximum and minimum degrees of 'a'.
    max_deg = sympy.degree(poly_in_a)
    # The `monoms()` method gives a list of tuples of the powers of the variables.
    powers_of_a = [m[0] for m in poly_in_a.monoms()]
    min_deg = min(powers_of_a)

    span_a = max_deg - min_deg

    print("\nFirst, we find the span of the polynomial:")
    print(f"The maximum power of 'a' is: {max_deg}")
    print(f"The minimum power of 'a' is: {min_deg}")
    print(f"The span is the difference: span_a(P) = {max_deg} - {min_deg} = {span_a}")

    # Step 3: Compute the lower bound using the inequality.
    lower_bound = (span_a / 2) + 1
    
    print("\nNow, we calculate the lower bound for s(9_23):")
    # Print the equation with the numbers substituted in.
    print(f"s(9_23) >= {span_a} / 2 + 1")
    print(f"s(9_23) >= {int(span_a / 2)} + 1")
    print(f"s(9_23) >= {int(lower_bound)}")

    print(f"\nThus, a lower bound for the minimum number of Seifert circles is {int(lower_bound)}.")

if __name__ == '__main__':
    solve_seifert_bound()