import sys
# The following code requires the 'spherogram' library.
# If you don't have it, please install it using: pip install spherogram
try:
    import spherogram
except ImportError:
    print("Error: The 'spherogram' library is required to run this code.")
    print("Please install it by running: pip install spherogram")
    sys.exit(1)

def solve_seifert_bound():
    """
    Calculates a lower bound for the minimum number of Seifert circles for a given knot
    using its HOMFLY polynomial.
    """
    knot_name = '9_23'
    try:
        knot = spherogram.Knot(knot_name)
    except KeyError:
        print(f"Error: Knot '{knot_name}' not found in the spherogram database.")
        return

    # 1. Compute the HOMFLY polynomial. The spherogram library uses variables l and m
    # where l corresponds to 'a' and m corresponds to 'z' in standard notation.
    homfly_poly = knot.homfly_polynomial()

    # 2. Extract the powers of the 'a' variable (which is the first variable, 'l').
    # The polynomial is stored as a dictionary where keys are tuples of exponents.
    exponents = homfly_poly.keys()
    if not exponents:
        print("Could not compute the polynomial or it is empty.")
        return
        
    a_powers = [exp[0] for exp in exponents]

    # 3. Find the minimum and maximum powers of 'a'.
    min_a_power = min(a_powers)
    max_a_power = max(a_powers)

    # 4. Calculate the span of the polynomial in 'a'.
    span_a = max_a_power - min_a_power

    # 5. Apply the Morton-Franks-Williams inequality to find the lower bound.
    # The formula is: s_min(K) >= span_a(P_K)/2 + 1
    lower_bound = (span_a / 2) + 1

    # 6. Print the detailed steps of the calculation.
    print(f"To find a lower bound for the minimum number of Seifert circles for the {knot_name} knot, we use the Morton-Franks-Williams inequality.")
    print("The inequality is: s_min(K) >= span_a(P(K))/2 + 1\n")

    print("Step 1: Compute the HOMFLY polynomial P(a, z).")
    # Note: Spherogram uses P(l,m). We will refer to it as P(a,z) for clarity.
    print(f"The HOMFLY polynomial for {knot_name} is: P(a, z) = {homfly_poly}\n")
    
    print("Step 2: Find the span of the polynomial in the variable 'a'.")
    print(f"The powers of 'a' in the polynomial are: {sorted(list(set(a_powers)))}")
    print(f"The minimum power of 'a' is: {min_a_power}")
    print(f"The maximum power of 'a' is: {max_a_power}")
    print(f"The span_a is max_power - min_power = {max_a_power} - {min_a_power} = {span_a}\n")

    print("Step 3: Calculate the lower bound using the inequality.")
    print(f"Lower Bound = (span_a / 2) + 1")
    print(f"Lower Bound = ({span_a} / 2) + 1")
    print(f"Lower Bound = {int(span_a / 2)} + 1")
    print(f"Lower Bound = {int(lower_bound)}\n")

    print(f"Thus, the lower bound for the minimum number of Seifert circles of the {knot_name} knot is {int(lower_bound)}.")

if __name__ == '__main__':
    solve_seifert_bound()