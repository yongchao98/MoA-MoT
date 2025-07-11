# Note: This script requires the 'spherogram' library.
# You can install it using pip: pip install spherogram
# Spherogram may have dependencies (like the SnapPy kernel) that also need to be installed.

import spherogram
import sympy

def solve_seifert_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles
    for the 9_23 knot using its HOMFLY polynomial.
    """
    try:
        # Define the knot
        knot = spherogram.Knot(9, 23)

        print(f"Finding a lower bound for the minimum number of Seifert circles of the {knot.name()} knot.")
        print("The method uses the Morton-Franks-Williams inequality on the HOMFLY polynomial.\n")

        # 1. Get the HOMFLY polynomial P(a, z)
        # The spherogram library computes this for us.
        homfly_poly = knot.homfly_polynomial()
        
        print("Step 1: Get the HOMFLY polynomial P(a, z).")
        print(f"For the {knot.name()} knot, the polynomial is:")
        print(f"P(a, z) = {homfly_poly}\n")
        
        # 2. Calculate the span of the variable 'a'
        # We can analyze the polynomial's terms to find the min and max powers of 'a'.
        # as_poly() with Laurent=True (default) handles negative exponents.
        poly_dict = homfly_poly.as_poly().as_dict()
        
        # The keys of the dictionary are tuples of exponents of the variables (a, z).
        # We extract the first element of each tuple, which is the power of 'a'.
        a_powers = [k[0] for k in poly_dict.keys()]
        
        min_a_power = min(a_powers)
        max_a_power = max(a_powers)
        span_a = max_a_power - min_a_power

        print("Step 2: Calculate the span of the 'a' variable in the polynomial.")
        print(f"The minimum power of 'a' found is: {min_a_power}")
        print(f"The maximum power of 'a' found is: {max_a_power}")
        print(f"The span is max_power - min_power = {max_a_power} - ({min_a_power}) = {span_a}\n")

        # 3. Apply the Morton-Franks-Williams inequality
        # s >= span_a(P) / 2 + 1
        # The span for a knot's HOMFLY polynomial is always even.
        lower_bound = (span_a // 2) + 1

        print("Step 3: Apply the Morton-Franks-Williams inequality: s >= span_a(P)/2 + 1")
        print("Substituting the calculated span:")
        print(f"s >= {span_a} / 2 + 1")
        print(f"s >= {span_a // 2} + 1")
        print(f"s >= {lower_bound}\n")
        
        print(f"The calculated lower bound for the minimum number of Seifert circles is {lower_bound}.")

    except ImportError:
        print("The 'spherogram' library is not installed.")
        print("Please install it by running: pip install spherogram")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("This may be due to a missing SnapPy kernel or other dependencies.")

if __name__ == "__main__":
    solve_seifert_bound()