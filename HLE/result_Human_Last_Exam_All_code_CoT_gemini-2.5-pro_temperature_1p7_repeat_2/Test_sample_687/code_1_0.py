import spherogram
import sympy

# This script requires the 'spherogram' library.
# You can install it via pip: pip install spherogram

# --- Step 1: Get the HOMFLY polynomial for the 9_23 knot ---
try:
    knot = spherogram.Knot('9_23')
    poly = knot.homfly_polynomial()
    print(f"The HOMFLY polynomial for the 9_23 knot is P(a,z) = {poly}")

    a, z = poly.variables

    # Expand the polynomial to easily access the powers of 'a'
    poly_expanded = sympy.expand(poly)
    
    # --- Step 2: Find the maximum and minimum powers of 'a' ---
    # Collect all powers of 'a' from the terms of the expanded polynomial
    coeffs = poly_expanded.as_coefficients_dict()
    a_powers = [term.as_powers_dict().get(a, 0) for term in coeffs.keys()]
    
    max_power = max(a_powers)
    min_power = min(a_powers)

    print(f"\nThe maximum power of 'a' in the polynomial is: {max_power}")
    print(f"The minimum power of 'a' in the polynomial is: {min_power}")

    # --- Step 3: Calculate the span of the polynomial in 'a' ---
    span_a = max_power - min_power
    print(f"The span of the polynomial is span_a = {max_power} - {min_power} = {span_a}")

    # --- Step 4: Calculate the lower bound for the number of Seifert circles ---
    # The formula is s(K) >= span_a / 2 + 1
    lower_bound = span_a / 2 + 1

    print("\nUsing the inequality s(K) >= span_a(P) / 2 + 1, we find the lower bound:")
    print(f"s(9_23) >= {span_a} / 2 + 1 = {int(lower_bound)}")
    print(f"\nThus, a lower bound for the minimum number of Seifert circles is {int(lower_bound)}.")

except ImportError:
    print("This script requires the 'spherogram' library.")
    print("Please install it using: pip install spherogram")
except Exception as e:
    print(f"An error occurred: {e}")
    print("Please ensure you have a working installation of spherogram and its dependencies.")
