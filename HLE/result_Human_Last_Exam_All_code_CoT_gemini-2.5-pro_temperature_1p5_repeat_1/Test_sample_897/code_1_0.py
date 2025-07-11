import spherogram
import sympy

def get_l_span(polynomial):
    """
    Computes the span of the 'l' variable in a sympy polynomial expression.
    This is the difference between the maximum and minimum degree of 'l'.
    """
    l, _ = polynomial.parent().gens()
    
    # Using as_coefficients_dict is robust for expressions with negative exponents
    term_dict = sympy.expand(polynomial).as_coefficients_dict()
    
    l_degrees = []
    for term, coeff in term_dict.items():
        if coeff != 0:
            powers = term.as_powers_dict()
            if l in powers:
                l_degrees.append(powers[l])
            else:
                # This term is constant with respect to l (e.g., a power of m)
                l_degrees.append(0)

    if not l_degrees:
        return 0
        
    return max(l_degrees) - min(l_degrees)

# --- Part 1: Analysis of K1 = 10_74 ---

print("--- Analyzing Knot K1 = 10_74 ---")
K1 = spherogram.Link('10_74')

# Calculate the HOMFLY polynomial.
# spherogram uses variables (l, m) where l corresponds to the variable 'a' or 'v' in standard notations.
poly_K1 = K1.homfly_polynomial()
print(f"The HOMFLY polynomial for K1 is: P(l,m) = {poly_K1}")

# Calculate the span of the polynomial in the 'l' variable
span_l_K1 = get_l_span(poly_K1)
print(f"The span of the polynomial in variable 'l' is: {span_l_K1}")

# Calculate the lower bound for the minimum number of Seifert circles
lower_bound_s_K1 = (span_l_K1 / 2) + 1
print(f"The lower bound for the min number of Seifert circles s(K1) is (span/2) + 1 = ({span_l_K1}/2) + 1 = {lower_bound_s_K1}")
print("-" * 30)


# --- Part 2: Analysis of K2 ---

print("--- Analyzing Knot K2 = closure((sigma_1^-1)^3 * sigma_2^-1) ---")
# This is a braid on 3 strands (since indices are 1 and 2).
# Braid word for (sigma_1^-1)^3 * sigma_2^-1 is [-1, -1, -1, -2]
K2 = spherogram.Link(braid_word=[-1, -1, -1, -2])

# The knot is represented by a 3-strand braid, so its braid index is at most 3.
braid_index_upper_bound_K2 = 3
print(f"K2 is given as a 3-braid, so its braid index b(K2) <= {braid_index_upper_bound_K2}.")
# The name of this knot in standard tables is 5_2.
print(f"(This knot is commonly known as {K2.name()})")

# Calculate the HOMFLY polynomial for K2
poly_K2 = K2.homfly_polynomial()
print(f"The HOMFLY polynomial for K2 is: P(l,m) = {poly_K2}")

# Calculate the span to find a lower bound for the braid index
span_l_K2 = get_l_span(poly_K2)
print(f"The span of the polynomial in variable 'l' is: {span_l_K2}")

# Calculate the lower bound for the braid index
lower_bound_b_K2 = (span_l_K2 / 2) + 1
print(f"The lower bound for the braid index b(K2) is (span/2) + 1 = ({span_l_K2}/2) + 1 = {lower_bound_b_K2}")

# Determine the exact braid index
print(f"Since b(K2) <= {braid_index_upper_bound_K2} and b(K2) >= {lower_bound_b_K2}, the braid index must be an integer.")
braid_index_K2 = int(lower_bound_b_K2)
print(f"Therefore, the braid index of K2 is exactly {braid_index_K2}.")
print("-" * 30)

# --- Part 3: Final Calculation ---
print("--- Final Calculation ---")
difference = braid_index_K2 - lower_bound_s_K1
print("The difference is (braid index of K2) - (lower bound for s(K1))")
print(f"= {braid_index_K2} - {lower_bound_s_K1} = {difference}")

<<< -2 >>>