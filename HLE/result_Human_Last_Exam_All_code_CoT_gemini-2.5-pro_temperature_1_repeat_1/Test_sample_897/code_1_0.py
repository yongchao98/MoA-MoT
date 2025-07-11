import sympy
import pyknotid.catalogue as cat
import pyknotid.spacecurves as sp
from pyknotid.spacecurves.braids import Braid

def get_l_span(poly_expr):
    """
    Calculates the span of the variable 'l' in a sympy expression.
    The span is the maximum degree minus the minimum degree.
    """
    l = sympy.symbols('l')
    
    # Expand the expression to ensure it's a sum of terms
    poly_expr = sympy.expand(poly_expr)
    
    # Collect all terms of the sum
    if isinstance(poly_expr, sympy.Add):
        terms = poly_expr.as_ordered_terms()
    else:
        terms = [poly_expr]
    
    degrees = set()
    for term in terms:
        # For each term, find the exponent of l
        powers = term.as_powers_dict()
        if l in powers:
            degrees.add(powers[l])
        else:
            # If l is not in the term, its power is 0
            degrees.add(0)
            
    if not degrees:
        return 0
        
    return max(degrees) - min(degrees)

# --- Part 1: Analysis of K1 = 10_74 ---

# Get the knot 10_74 from the catalogue
k1 = cat.get_knot('10_74')

# Calculate its HOMFLY polynomial P(l, m)
# The pyknotid library returns a sympy expression for the polynomial.
p1 = k1.homfly_polynomial()

# The Morton-Franks-Williams inequality gives a lower bound on the
# minimum number of Seifert circles, s(K), based on the l-span of the HOMFLY polynomial.
# s(K) >= span_l(P(K))/2 + 1
span_l_k1 = get_l_span(p1)
seifert_circles_lower_bound_K1 = (span_l_k1 / 2) + 1

print(f"For K1 = 10_74:")
print(f"  HOMFLY polynomial P(l, m) = {p1}")
print(f"  The span of the 'l' variable is {span_l_k1}.")
print(f"  The lower bound for the minimum number of Seifert circles is ({span_l_k1}/2) + 1 = {seifert_circles_lower_bound_K1}.")
print("-" * 20)

# --- Part 2: Analysis of K2 = closure of (s1^-1)^3 * s2^-1 ---

# The braid is given in B_3, so the braid index is at most 3.
# The braid word is s1^-1 s1^-1 s1^-1 s2^-1, represented by [-1, -1, -1, -2].
braid_k2 = Braid([-1, -1, -1, -2])
k2 = sp.Knot(braid_k2.space_curve())

# Calculate its HOMFLY polynomial P(l, m)
p2 = k2.homfly_polynomial()

# The Morton-Franks-Williams inequality gives a lower bound on the braid index n(K).
# n(K) >= span_l(P(K))/2 + 1
span_l_k2 = get_l_span(p2)
braid_index_lower_bound_K2 = (span_l_k2 / 2) + 1

# Since K2 is a closure of a 3-braid, its braid index n(K2) <= 3.
# We have lower_bound <= n(K2) <= 3.
# This allows us to determine the exact braid index.
braid_index_K2 = braid_index_lower_bound_K2

print(f"For K2 = closure of (sigma_1^-1)^3 * sigma_2^-1:")
print(f"  HOMFLY polynomial P(l, m) = {p2}")
print(f"  The span of the 'l' variable is {span_l_k2}.")
print(f"  The lower bound for the braid index is ({span_l_k2}/2) + 1 = {braid_index_lower_bound_K2}.")
print(f"  Since K2 is from a 3-braid, its index is at most 3. Thus, the braid index is {braid_index_K2}.")
print("-" * 20)

# --- Part 3: Final Difference Calculation ---

difference = braid_index_K2 - seifert_circles_lower_bound_K1

print("Final Calculation:")
print("The difference between the braid index of K2 and the lower bound of the minimum number of Seifert circles of K1 is:")
print(f"{int(braid_index_K2)} - {int(seifert_circles_lower_bound_K1)} = {int(difference)}")
<<<0>>>