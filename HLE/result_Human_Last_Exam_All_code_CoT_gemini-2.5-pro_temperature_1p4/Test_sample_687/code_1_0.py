# The solution requires the 'spherogram' and 'sympy' libraries.
# You can install them using pip:
# pip install spherogram
# pip install sympy

import spherogram
import sympy

# --- Introduction ---
print("This script calculates a lower bound for the minimum number of Seifert circles, s(K), for the 9_23 knot.")
print("We use the inequality derived from the HOMFLY polynomial P(a,z): s(K) >= span_z(P)/2 + 1.\n")

# --- Step 1: Get the knot and its HOMFLY polynomial ---
print("Step 1: Compute the HOMFLY polynomial for the 9_23 knot.")
try:
    # In spherogram, knots are identified by their Alexander-Briggs notation, e.g., (9, 23) for 9_23.
    knot = spherogram.Knot(9, 23)

    # Compute the HOMFLY polynomial. Spherogram uses variables L and M.
    # The span of the second variable, M, is what we need for the inequality.
    poly = knot.homfly_polynomial()
    print(f"Successfully computed the HOMFLY polynomial for {knot}.")

except ImportError:
    print("Error: The 'spherogram' library is required but not found.")
    print("Please install it using: pip install spherogram")
    exit()
except Exception as e:
    print(f"An error occurred: {e}")
    exit()

# --- Step 2: Calculate the span of the polynomial ---
print("\nStep 2: Determine the span of the polynomial in its second variable (M).")

# The poly object from spherogram allows access to its terms.
# Each term is a tuple of the form ((exponent_L, exponent_M), coefficient).
all_terms = poly.terms()
m_exponents = [term[0][1] for term in all_terms]

min_deg_m = min(m_exponents)
max_deg_m = max(m_exponents)
span_m = max_deg_m - min_deg_m

print(f"The powers of the second variable found in the polynomial are: {sorted(list(set(m_exponents)))}")
print(f"The minimum power is: {min_deg_m}")
print(f"The maximum power is: {max_deg_m}")
print(f"The span is the difference: span = {max_deg_m} - {min_deg_m} = {span_m}")

# --- Step 3: Apply the inequality to find the lower bound ---
print("\nStep 3: Apply the inequality to find the lower bound for the number of Seifert circles.")

# The formula is s(K) >= span / 2 + 1. The result must be an integer.
s_lower_bound = (span_m // 2) + 1

print(f"The formula is: s(9_23) >= span / 2 + 1")
print(f"Substituting the span value: s(9_23) >= {span_m} / 2 + 1 = {span_m // 2} + 1 = {s_lower_bound}")
print("\n----------------------------------------------------------------")
print(f"A lower bound for the minimum number of Seifert circles of the 9_23 knot is {s_lower_bound}.")
print("----------------------------------------------------------------")
