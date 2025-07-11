# This script requires the 'pyknotid' and 'sympy' libraries.
# You can install them using: pip install pyknotid sympy

import sympy
from pyknotid.catalogue import get_knot

try:
    # Step 1: Identify the knot, 9_23.
    knot = get_knot(9, 23)
    print(f"Finding the lower bound for the minimum number of Seifert circles for the {knot.identifier} knot.")

    # Step 2: Calculate its HOMFLY polynomial.
    # The pyknotid library returns a sympy expression in variables 'l' and 'm'.
    homfly_poly = knot.homfly_polynomial()
    print(f"\nThe HOMFLY polynomial P(l, m) is: {homfly_poly}")

    # The bound for the Seifert circle number s(K) is derived from the span of the polynomial.
    # The inequality is s(K) >= span_a(P)/2 + 1.
    # The span in the variable 'l' of P(l,m) is equal to the span in 'a' of P(a,z).

    # Step 3: Find the span of the polynomial in the variable 'l'.
    # We first need to identify the variables and then get all the exponents for 'l'.
    l, m = sympy.symbols('l, m')
    poly_in_l_m = sympy.poly(homfly_poly, l, m)
    
    # Extract all exponents of the 'l' variable from the polynomial terms.
    l_exponents = [term[0][0] for term in poly_in_l_m.terms()]

    max_l_degree = max(l_exponents)
    min_l_degree = min(l_exponents)
    span = max_l_degree - min_l_degree

    print(f"\nThe highest power of 'l' is: {max_l_degree}")
    print(f"The lowest power of 'l' is: {min_l_degree}")
    
    # Step 4: Calculate the lower bound using the formula.
    # We print each number in the final equation.
    print("\nThe span of the polynomial is calculated as max_degree - min_degree.")
    print(f"span = {max_l_degree} - {min_l_degree} = {span}")

    lower_bound = span / 2 + 1
    
    print("\nThe lower bound for the number of Seifert circles is span/2 + 1.")
    print(f"Lower Bound = {span}/2 + 1 = {int(span/2)} + 1 = {int(lower_bound)}")

    print(f"\nConclusion: The minimum number of Seifert circles for the 9_23 knot is at least {int(lower_bound)}.")

except ImportError:
    print("This script requires 'pyknotid' and 'sympy'. Please install them by running: pip install pyknotid sympy")

except Exception as e:
    print(f"An error occurred: {e}")
    print("Please ensure you have an internet connection for the knot catalogue, or that the libraries are installed correctly.")
