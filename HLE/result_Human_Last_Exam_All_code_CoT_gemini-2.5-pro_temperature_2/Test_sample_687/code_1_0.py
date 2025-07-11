import spherogram
import sympy

def solve_knot_problem():
    """
    This function calculates a lower bound for the minimum number of Seifert
    circles of the 9_23 knot using its HOMFLY polynomial.
    """
    # 1. Identify the knot
    knot_name = '9_23'
    try:
        k = spherogram.Knot(knot_name)
    except ImportError:
        print("Spherogram library not found. Please install it using 'pip install spherogram'.")
        return
    except Exception as e:
        print(f"An error occurred while loading the knot {knot_name}: {e}")
        return

    # 2. Compute the HOMFLY polynomial
    # The spherogram library uses variables l and m. The span of l corresponds
    # to the span of v in other conventions.
    poly = k.homfly_polynomial()
    l_var = poly.gens[0]

    print(f"The knot is {knot_name}.")
    # print(f"The HOMFLY polynomial is P(l, m) = {poly}")

    # 3. Calculate the span of the 'l' variable
    try:
        # Treat the expression as a polynomial in the variable 'l'
        poly_in_l = sympy.Poly(poly, l_var)
        
        # Find the maximum and minimum degrees of l
        max_deg = sympy.degree(poly_in_l)
        min_deg = sympy.lowest_degree(poly_in_l)
        
        span = max_deg - min_deg

    except Exception as e:
        print(f"Could not compute the degrees from the polynomial. Error: {e}")
        return

    # 4. Compute the lower bound for the number of Seifert circles
    # The formula is: s_min(K) >= span/2 + 1
    lower_bound = span / 2 + 1

    print("To find the lower bound for the minimum number of Seifert circles, we use the formula:")
    print("s_min >= span/2 + 1")
    print("\nFirst, we find the span of the HOMFLY polynomial's first variable (l).")
    print(f"The maximum degree of l is: {max_deg}")
    print(f"The minimum degree of l is: {min_deg}")
    print(f"The span of l is: {max_deg} - ({min_deg}) = {span}")

    print("\nNow we can calculate the lower bound:")
    print(f"Lower Bound = {span} / 2 + 1 = {int(lower_bound)}")
    print(f"\nThus, a lower bound for the minimum number of Seifert circles of the {knot_name} knot is {int(lower_bound)}.")

solve_knot_problem()