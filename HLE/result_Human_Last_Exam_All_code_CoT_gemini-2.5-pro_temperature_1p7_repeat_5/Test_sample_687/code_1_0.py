import spherogram
import sympy

def solve_seifert_circle_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles of a knot
    using its HOMFLY polynomial.
    """
    # The knot we are interested in is 9_23.
    knot_name = '9_23'
    try:
        knot = spherogram.Knot(knot_name)
    except Exception as e:
        print(f"Could not create knot {knot_name}. Make sure the 'spherogram' library is installed (`pip install spherogram`).")
        print(f"Error: {e}")
        return

    # The Morton-Franks-Williams inequality provides the lower bound for the
    # minimum number of Seifert circles, s(K), from the HOMFLY polynomial P(a, z):
    # s(K) >= span_z(P) / 2 + 1

    # 1. Compute the HOMFLY polynomial for the 9_23 knot.
    # The spherogram library computes the HOMFLY-PT polynomial in variables (v, z).
    # The span in 'z' is invariant across common HOMFLY polynomial conventions.
    homfly_poly = knot.homfly_polynomial()

    # 2. Find the span of the polynomial in the variable 'z'.
    z = sympy.Symbol('z')
    
    # Treat the expression as a polynomial in 'z' to find the degrees.
    poly_in_z = sympy.Poly(homfly_poly, z)
    
    # Get the highest and lowest powers of 'z' from the polynomial.
    max_degree_z = sympy.degree(poly_in_z, gen=z)
    min_degree_z = poly_in_z.min_degree()

    # The span is the difference between the maximum and minimum degrees.
    span_z = max_degree_z - min_degree_z

    # 3. Apply the inequality to find the lower bound.
    lower_bound = span_z / 2 + 1

    # 4. Print the final result, showing the calculation step-by-step.
    print(f"To find a lower bound for the minimum number of Seifert circles, s(K), for the {knot_name} knot, we use the inequality:")
    print("s(K) >= span_z(P) / 2 + 1")
    print(f"The maximum degree of z in the HOMFLY polynomial is {max_degree_z}.")
    print(f"The minimum degree of z is {min_degree_z}.")
    print(f"The span in z is the difference: {max_degree_z} - {min_degree_z} = {span_z}.")
    print("Plugging the span into the formula, we get the lower bound:")
    # The final equation with all numbers is printed below
    print(f"s(K) >= {span_z} / 2 + 1 = {int(lower_bound)}")
    print(f"\nA lower bound for the minimum number of Seifert circles of the {knot_name} knot is {int(lower_bound)}.")

# Run the solver
solve_seifert_circle_bound()