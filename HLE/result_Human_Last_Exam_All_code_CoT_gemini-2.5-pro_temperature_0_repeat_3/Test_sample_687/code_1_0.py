import spherogram
import sympy

def solve_knot_problem():
    """
    Calculates a lower bound for the minimum number of Seifert circles
    for the 9_23 knot using its HOMFLY polynomial.
    """
    try:
        # 1. Create a Knot object for 9_23
        knot = spherogram.Knot('9_23')

        # 2. Calculate its HOMFLY-PT polynomial
        poly = knot.homfly_pt_polynomial()
        
        print(f"The HOMFLY polynomial P(a, z) for the 9_23 knot is: {poly}")

        # The variables are 'a' and 'z'
        a, z = sympy.symbols('a, z')

        # 3. Convert the expression to a polynomial in z to easily get the degrees
        poly_in_z = sympy.Poly(poly, z)

        # Get all the degrees of z present in the polynomial
        # monoms() returns a list of tuples, where each tuple contains the exponent of z
        z_degrees = [m[0] for m in poly_in_z.monoms()]

        # 4. Find the maximum and minimum degrees of z
        max_deg = max(z_degrees)
        min_deg = min(z_degrees)

        print(f"\nThe degrees of the variable z in the polynomial are: {sorted(z_degrees)}")
        print(f"The maximum degree of z is {max_deg}.")
        print(f"The minimum degree of z is {min_deg}.")

        # 5. Calculate the span of the polynomial in z
        span = max_deg - min_deg
        print(f"The span of the polynomial in z is calculated as: {max_deg} - {min_deg} = {span}")

        # 6. Use the inequality s(K) >= span_z/2 + 1 to find the lower bound
        lower_bound = (span / 2) + 1
        
        print("\nA lower bound for the minimum number of Seifert circles, s(K), is given by the inequality:")
        print("s(K) >= span_z / 2 + 1")
        print("\nSubstituting the calculated span:")
        # The final equation with each number printed
        print(f"s(K) >= {span} / 2 + 1 = {int(span/2)} + 1 = {int(lower_bound)}")
        
        print(f"\nTherefore, a lower bound for the minimum number of Seifert circles of the 9_23 knot is {int(lower_bound)}.")

    except ImportError:
        print("Please install the 'spherogram' and 'sympy' libraries to run this code.")
        print("You can install them using: pip install spherogram sympy")
    except Exception as e:
        print(f"An error occurred: {e}")

solve_knot_problem()