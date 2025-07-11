import re

def get_polynomial_span(poly_str, var):
    """
    A simple function to find the span of a polynomial string.
    It finds all occurrences of the variable with an exponent
    and computes the difference between the max and min exponent.
    Example: get_polynomial_span("x**-2 + x**-4", "x") returns 2.
    """
    pattern = re.compile(r"\%s\*\*([\-0-9]+)" % var)
    exponents = [int(p) for p in pattern.findall(poly_str)]
    if not exponents:
        return 0
    return max(exponents) - min(exponents)

# Part 1: Calculation for K1 = 10_74 knot

print("--- Part 1: Lower bound of Seifert circles for K1 (10_74) ---")

# The HOMFLY polynomial for the 10_74 knot is
# P(a, z) = a^{-6} - a^{-6}z^2 + 2a^{-4}z^2 - a^{-4}z^4 + a^{-2}z^4
# We can rewrite this grouped by powers of 'a':
# P(a, z) = (1-z^2)*a^{-6} + (2z^2-z^4)*a^{-4} + (z^4)*a^{-2}
homfly_10_74 = "(1-z**2)*a**-6 + (2*z**2-z**4)*a**-4 + (z**4)*a**-2"
print(f"The HOMFLY polynomial for 10_74 is: {homfly_10_74}")

# The powers of 'a' in the polynomial are -6, -4, and -2.
max_deg_a = -2
min_deg_a = -6

# Calculate the span of the polynomial in variable 'a'.
span_a = max_deg_a - min_deg_a
print(f"The maximum degree of 'a' is {max_deg_a}.")
print(f"The minimum degree of 'a' is {min_deg_a}.")
print(f"The span of the HOMFLY polynomial in 'a' is: {max_deg_a} - ({min_deg_a}) = {span_a}")

# The Morton-Franks-Williams inequality gives a lower bound on the
# minimum number of Seifert circles, s(K).
# s(K) >= span_a(P)/2 + 1
seifert_circles_lower_bound = span_a / 2 + 1
print(f"The lower bound for the minimum number of Seifert circles is: {span_a}/2 + 1 = {int(seifert_circles_lower_bound)}")

# Part 2: Calculation for K2 = closure of (sigma_1^-1)^3 * sigma_2^-1

print("\n--- Part 2: Braid index for K2 ---")
# The knot K2 is the closure of the braid (sigma_1^-1)^3 * sigma_2^-1
# from the 3-strand braid group B_3. This knot is identified as 5_2.
print("K2 is the closure of a 3-strand braid, identified as the knot 5_2.")
print("This implies its braid index b(K2) must be less than or equal to 3.")

# We use the Jones polynomial of 5_2 to find a lower bound for its braid index.
# The Jones polynomial for 5_2 is V(t) = t^-2 + t^-3 - t^-4.
jones_5_2 = "t**-2 + t**-3 - t**-4"
print(f"The Jones polynomial for 5_2 is: {jones_5_2}")

# The powers of 't' in the polynomial are -2, -3, and -4.
max_deg_t = -2
min_deg_t = -4
span_t = max_deg_t - min_deg_t
print(f"The maximum degree of 't' is {max_deg_t}.")
print(f"The minimum degree of 't' is {min_deg_t}.")
print(f"The span of the Jones polynomial in 't' is: {max_deg_t} - ({min_deg_t}) = {span_t}")

# The Morton-Franks-Williams inequality for the Jones polynomial is b(K) >= span_t(V) + 1.
braid_index_lower_bound = span_t + 1
print(f"The lower bound for the braid index is: {span_t} + 1 = {braid_index_lower_bound}")

# Since b(K2) <= 3 and b(K2) >= 3, the braid index is exactly 3.
braid_index_K2 = 3
print(f"Combining the bounds, the braid index of K2 is exactly {braid_index_K2}.")

# Part 3: Final Difference
print("\n--- Part 3: Final Calculation ---")
difference = braid_index_K2 - seifert_circles_lower_bound
print("The question asks for the difference between the braid index of K2 and the lower bound of Seifert circles of K1.")
print(f"Result = (Braid index of K2) - (Lower bound for Seifert circles of K1)")
print(f"Result = {braid_index_K2} - {int(seifert_circles_lower_bound)} = {int(difference)}")
