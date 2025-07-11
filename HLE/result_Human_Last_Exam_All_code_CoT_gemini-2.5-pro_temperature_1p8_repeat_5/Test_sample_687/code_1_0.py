import spherogram

def solve():
    """
    Calculates a lower bound for the minimum number of Seifert circles
    of the 9_23 knot using its HOMFLY polynomial.
    """
    try:
        # 1. Define the knot
        knot_name = '9_23'
        k = spherogram.Knot(knot_name)

        # 2. Calculate the HOMFLY polynomial P(v, z)
        # The spherogram library calculates the polynomial P(v,z) satisfying
        # the skein relation v^{-1}P(L_+) - vP(L_-) = zP(L_0).
        homfly_poly = k.homfly_polynomial()

        # 3. Find the span of the polynomial in the variable 'v'
        # The polynomial is a SageMath symbolic expression. We extract the exponents of v.
        # The variables are ordered as (v, z), so the first element of each exponent tuple is for 'v'.
        v_exponents = [exp[0] for exp in homfly_poly.exponents()]
        
        if not v_exponents:
            print(f"Could not determine exponents for the HOMFLY polynomial of {knot_name}.")
            return

        min_v_degree = min(v_exponents)
        max_v_degree = max(v_exponents)
        span_v = max_v_degree - min_v_degree

        # 4. Calculate the lower bound for the number of Seifert circles, s(K)
        # The bound is given by the formula: s(K) >= span_v(P)/2 + 1
        # The span for a knot is always an even integer, so division by 2 is exact.
        seifert_circles_lower_bound = (span_v / 2) + 1

        # 5. Print the results clearly, showing each number in the final equation.
        print(f"The knot is {knot_name}.")
        print(f"Its HOMFLY polynomial P(v, z) is: {homfly_poly}")
        print(f"The minimum degree of 'v' in the polynomial is: {min_v_degree}")
        print(f"The maximum degree of 'v' in the polynomial is: {max_v_degree}")
        print(f"The span of the HOMFLY polynomial in 'v' is max_degree - min_degree.")
        print(f"span_v(P) = {max_v_degree} - ({min_v_degree}) = {span_v}")
        print("\nA lower bound for the minimum number of Seifert circles, s(K), is derived from the Morton-Franks-Williams inequality:")
        print("s(K) >= span_v(P)/2 + 1")
        print(f"s(K) >= {span_v}/2 + 1")
        print(f"s(K) >= {span_v/2} + 1")
        print(f"s(K) >= {seifert_circles_lower_bound}")
        print(f"\nThus, a lower bound for the minimum number of Seifert circles of the {knot_name} knot is {int(seifert_circles_lower_bound)}.")

    except ImportError:
        print("This script requires the 'spherogram' library.")
        print("Please install it, e.g., using 'pip install spherogram'.")
    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure that spherogram and its dependencies (like SageMath) are correctly installed.")

if __name__ == "__main__":
    solve()