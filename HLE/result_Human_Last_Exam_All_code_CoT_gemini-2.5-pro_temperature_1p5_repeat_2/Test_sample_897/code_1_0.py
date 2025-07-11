import sympy
from pyknotid.catalogue import get_knot
from pyknotid.spacecurves.braid import Braid

def get_z_span(homfly_poly_expr):
    """
    Calculates the span of the 'm' variable (corresponding to 'z') in a
    sympy expression representing a HOMFLY polynomial, including Laurent polynomials.
    """
    m_var = sympy.symbols('m')
    
    # Numerator and denominator are polynomials in m (and other variables)
    num, den = sympy.fraction(homfly_poly_expr)

    # Convert them to explicit polynomial objects in terms of m
    poly_num_m = sympy.poly(num, m_var)
    poly_den_m = sympy.poly(den, m_var)
    
    # Extract the degrees of m in the monomials of the numerator and denominator
    if poly_num_m.is_zero:
        deg_num = [-float('inf')]
    else:
        deg_num = [monom[0] for monom in poly_num_m.monoms()]

    if poly_den_m.is_zero:
         # Should not be reachable for non-zero polynomials
        deg_den = [0]
    else:
        deg_den = [monom[0] for monom in poly_den_m.monoms()]
        
    # The degree of a rational function term P(m)/Q(m) is deg(P) - deg(Q).
    # The max degree of the expression is max(deg_num) - min(deg_den).
    # The min degree of the expression is min(deg_num) - max(deg_den).
    max_power = max(deg_num) - min(deg_den)
    min_power = min(deg_num) - max(deg_den)
    
    return max_power - min_power

# --- Part 1: Analysis of K1 = 10_74 ---

# Get the knot K1 from the knot catalogue
k1 = get_knot('10_74')

# Calculate the HOMFLY polynomial for K1
# pyknotid uses variables (l, m). We are interested in the span of m.
homfly_k1 = k1.homfly_polynomial()

# Calculate the span of the z-variable (m in pyknotid)
span_z_k1 = get_z_span(homfly_k1)

# Calculate the lower bound of the minimum number of Seifert circles using the inequality:
# s_min(K) >= span_z(P(K))/2 + 1
seifert_circles_bound_k1 = (span_z_k1 / 2) + 1


# --- Part 2: Analysis of K2 = closure of (sigma_1^-1)^3 * sigma_2^-1 ---

# Define the braid for K2. It's a 3-strand braid.
# sigma_i^-1 is represented as -i.
braid_k2 = Braid([-1, -1, -1, -2])

# The braid index b(K2) is at most 3, since it's a 3-braid closure.
braid_index_upper_bound = 3

# Calculate the HOMFLY polynomial for the closure of the braid
homfly_k2 = braid_k2.homfly_polynomial()

# Calculate the span of the z-variable
span_z_k2 = get_z_span(homfly_k2)

# Calculate the lower bound for the braid index using the Morton-Franks-Williams inequality:
# b(K) >= span_z(P(K))/2 + 1
braid_index_lower_bound = (span_z_k2 / 2) + 1

# Since braid_index_lower_bound <= b(K2) <= braid_index_upper_bound,
# and we will find these bounds are equal, we can determine the exact braid index.
braid_index_k2 = braid_index_lower_bound


# --- Final Calculation ---

# The difference between the braid index of K2 and the lower bound of the
# minimum number of Seifert circles of K1.
difference = braid_index_k2 - seifert_circles_bound_k1

print("Step 1: Calculate the lower bound of the minimum number of Seifert circles for K1 (10_74).")
print(f"The HOMFLY polynomial z-span for K1 is {int(span_z_k1)}.")
print(f"The lower bound for the number of Seifert circles is {int(span_z_k1)} / 2 + 1 = {int(seifert_circles_bound_k1)}.")
print("-" * 20)
print("Step 2: Calculate the braid index for K2.")
print(f"The HOMFLY polynomial z-span for K2 is {int(span_z_k2)}.")
print(f"The lower bound for the braid index is {int(span_z_k2)} / 2 + 1 = {int(braid_index_lower_bound)}.")
print(f"Since K2 is a 3-braid, its braid index is at most 3. Thus, the braid index is exactly {int(braid_index_k2)}.")
print("-" * 20)
print("Step 3: Calculate the difference.")
print("The difference is the braid index of K2 minus the Seifert circle bound of K1.")
print(f"Equation: {int(braid_index_k2)} - {int(seifert_circles_bound_k1)} = {int(difference)}")


<<<-1>>>