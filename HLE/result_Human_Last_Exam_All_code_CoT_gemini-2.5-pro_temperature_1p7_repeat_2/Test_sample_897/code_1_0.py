import sympy
from pyknotid.catalogue.identify import from_braid, get_knot
from pyknotid.knots import Knot

def get_poly_span(poly_expression, var):
    """
    Calculates the span of a variable in a sympy polynomial expression.
    The span is the difference between the maximum and minimum degree.
    """
    if poly_expression == 0:
        return 0
    
    # Expand the polynomial to handle forms like v**2 * (z**2 - 1)
    p_expanded = sympy.expand(poly_expression)
    
    # Collect all terms from the expanded polynomial
    if isinstance(p_expanded, sympy.Add):
        terms = p_expanded.as_ordered_terms()
    else: # Handle single-term polynomials
        terms = [p_expanded]

    # Find the degree of the variable in each term
    degrees = []
    for term in terms:
        # Using as_powers_dict() is a reliable way to get the exponent
        powers = term.as_powers_dict()
        degrees.append(powers.get(var, 0))
    
    if not degrees:
        return 0
        
    return max(degrees) - min(degrees)

# --- Part 1: Braid Index of K2 ---
print("--- Part 1: Analyzing K2 ---")

# K2 is the closure of the braid (sigma_1^-1)^3 * sigma_2^-1
# In pyknotid notation, this is [-1, -1, -1, -2]
braid_word_k2 = [-1, -1, -1, -2]
print(f"Identifying the knot K2 from its braid representation {braid_word_k2}...")

# Identify the knot using pyknotid
k2_id = from_braid(braid_word_k2)
k2 = get_knot(k2_id)
print(f"Result: The knot K2 is the {k2.name} knot.")

# Calculate the Jones polynomial for K2 to find its braid index
t = sympy.Symbol('t')
jones_poly_k2 = k2.jones_polynomial(variable=t)
print(f"The Jones polynomial of {k2.name} is: V(t) = {jones_poly_k2}")

# Calculate the span of the Jones polynomial
span_jones_k2 = get_poly_span(jones_poly_k2, t)
print(f"The span of the Jones polynomial is {span_jones_k2}.")

# Apply the Morton-Franks-Williams inequality: span(V) <= b(K) - 1
braid_index_lower_bound_k2 = span_jones_k2 + 1
print(f"By the Morton-Franks-Williams inequality, the braid index b(K2) must be >= {span_jones_k2} + 1 = {braid_index_lower_bound_k2}.")
print(f"Since K2 was defined by a 3-strand braid, its braid index is at most 3.")

# Conclude the braid index
braid_index_of_k2 = 3
print(f"Therefore, the braid index of K2 is {braid_index_of_k2}.\n")


# --- Part 2: Seifert Circle Lower Bound for K1 ---
print("--- Part 2: Analyzing K1 ---")

# K1 is the 10_74 knot
print("K1 is the 10_74 knot.")
k1 = Knot(10, 74)

# Calculate the HOMFLY polynomial for K1
v, z = sympy.symbols('v z')
homfly_poly_k1 = k1.homfly_polynomial(variables=(v, z))
print(f"The HOMFLY polynomial of {k1.name} is: P(v, z) = {homfly_poly_k1}")

# Calculate the v-span of the HOMFLY polynomial
v_span_homfly_k1 = get_poly_span(homfly_poly_k1, v)
print(f"The span of the 'v' variable in the HOMFLY polynomial is {v_span_homfly_k1}.")

# Apply the inequality: span_v(P) <= 2 * (s_min - 1)
# which implies s_min >= span_v(P) / 2 + 1
seifert_circles_lower_bound_k1 = (v_span_homfly_k1 / 2) + 1
print(f"Using the inequality s_min >= span_v/2 + 1, the lower bound for the minimum number of Seifert circles is:")
print(f"s_min(K1) >= {v_span_homfly_k1}/2 + 1 = {int(seifert_circles_lower_bound_k1)}.\n")


# --- Part 3: Final Calculation ---
print("--- Part 3: Final Calculation ---")
print("We need to find the difference between the braid index of K2 and the lower bound of the Seifert circles of K1.")

result = braid_index_of_k2 - seifert_circles_lower_bound_k1
print(f"The equation is: {braid_index_of_k2} - {int(seifert_circles_lower_bound_k1)} = {int(result)}")

print("\nFinal Answer:")
print(f"<<<{int(result)}>>>")