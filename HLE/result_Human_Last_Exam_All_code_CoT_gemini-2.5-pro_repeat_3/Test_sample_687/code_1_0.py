# First, ensure you have the necessary libraries installed.
# You can install them by running: pip install snek-knot sympy

import snek_knot
import sympy

def solve_knot_problem():
    """
    Calculates a lower bound for the minimum number of Seifert circles
    of the 9_23 knot using its HOMFLY polynomial.
    """
    try:
        # Step 1: Get the knot data for 9_23
        knot_name = "9_23"
        knot = snek_knot.Knot(knot_name)

        # The HOMFLY polynomial is often denoted P(a, z).
        # The 'snek-knot' library uses variables L and M.
        # L corresponds to 'a' and M corresponds to 'z' in the standard P(a, z) notation.
        # The library returns the polynomial as a dictionary where keys are the powers
        # of 'a' (as integers) and values are the coefficient polynomials in 'z' (as sympy expressions).
        homfly_poly_dict = knot.homfly_pt_polynomial()

        # Step 2: Identify the powers of the variable 'a'
        a_powers = list(homfly_poly_dict.keys())

        # Step 3: Find the maximum and minimum powers of 'a'
        max_power_a = max(a_powers)
        min_power_a = min(a_powers)

        # Step 4: Calculate the span of 'a'
        # span_a = max_degree(a) - min_degree(a)
        span_a = max_power_a - min_power_a

        # Step 5: Apply the theorem to find the lower bound for the number of Seifert circles
        # s(K) >= span_a(P) + 1
        lower_bound = span_a + 1
        
        # --- Output the results step-by-step ---
        print(f"Finding a lower bound for the minimum number of Seifert circles of the {knot_name} knot.")
        print("-" * 70)
        
        print("The HOMFLY polynomial P(a, z) for this knot has the following structure:")
        # We replace the library's variable M with z for clarity
        z = sympy.Symbol('z')
        for power in sorted(homfly_poly_dict.keys(), reverse=True):
            coeff_poly = homfly_poly_dict[power].subs(sympy.Symbol('M'), z)
            print(f"  Term with a^{power}: ({coeff_poly})")
        
        print("\nStep 1: Find the maximum and minimum powers of the variable 'a'.")
        print(f"The powers of 'a' present in the polynomial are: {sorted(a_powers)}")
        print(f"The maximum power of 'a' is: {max_power_a}")
        print(f"The minimum power of 'a' is: {min_power_a}")

        print("\nStep 2: Calculate the span of the polynomial in 'a'.")
        print(f"The span is the difference between the maximum and minimum powers.")
        # The final equation with each number explicitly shown
        print(f"span_a = {max_power_a} - ({min_power_a}) = {span_a}")

        print("\nStep 3: Apply the theorem s(K) >= span_a + 1 to find the lower bound.")
        # The final equation with each number explicitly shown
        print(f"Lower Bound = {span_a} + 1 = {lower_bound}")
        
        print("-" * 70)
        print(f"The calculated lower bound for the minimum number of Seifert circles is {lower_bound}.")

    except ImportError:
        print("Error: The 'snek-knot' or 'sympy' library is not installed.")
        print("Please install them by running: pip install snek-knot sympy")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    solve_knot_problem()