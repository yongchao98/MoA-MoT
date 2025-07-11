import spherogram
import sympy

# This script requires the 'spherogram' and 'sympy' libraries.
# You can install them using: pip install spherogram sympy

# Part 1: Analyze K2
# K2 is the closure of the braid (sigma_1^{-1})^3 * sigma_2^{-1} in the 3-strand braid group B_3.
# In the notation used by spherogram, this braid is represented by the list [-1, -1, -1, -2].
braid_representation_K2 = [-1, -1, -1, -2]
K2 = spherogram.Link(braid_representation_K2)

# The braid_index() method in spherogram computes the braid index of the knot.
braid_index_K2 = K2.braid_index()
print(f"The knot K2 is the closure of the given braid. It is identified as the knot {K2.identify()[0]}.")
print(f"The braid index of K2 is {braid_index_K2}.")

# Part 2: Analyze K1
# K1 is the knot 10_74.
K1 = spherogram.Link('10_74')
print(f"\nThe knot K1 is {K1.name}.")

# Compute the HOMFLY polynomial for K1. Spherogram uses variables (l, m).
# The second variable 'm' corresponds to the 'z' variable in the standard HOMFLY polynomial P(a, z)
# for the purpose of calculating the z-span.
homfly_K1 = K1.homfly_polynomial()

# We can display the polynomial for clarity using sympy.
l, m = sympy.symbols('l, m')
homfly_K1_sympy = homfly_K1.to_sympy(l, m)
print(f"The HOMFLY polynomial of K1 is P(l,m) = {homfly_K1_sympy}")


# To find the lower bound for the number of Seifert circles, we need the span of the polynomial
# in its second variable ('m' in spherogram, 'z' in standard notation).
# The span is the difference between the maximum and minimum exponents.
m_exponents = [e[1] for e in homfly_K1.exponents()]
max_m_exp = max(m_exponents)
min_m_exp = min(m_exponents)
z_span = max_m_exp - min_m_exp

print(f"The maximum power of the second variable 'm' is {max_m_exp}.")
print(f"The minimum power of the second variable 'm' is {min_m_exp}.")
print(f"The span of the polynomial in the second variable is {max_m_exp} - ({min_m_exp}) = {z_span}.")

# The lower bound for the minimum number of Seifert circles is span_z + 1.
seifert_circles_lower_bound_K1 = z_span + 1
print(f"The lower bound for the minimum number of Seifert circles of K1 is {z_span} + 1 = {seifert_circles_lower_bound_K1}.")

# Part 3: Calculate the difference
# The question asks for the difference between the braid index of K2 and the lower bound for K1.
difference = braid_index_K2 - seifert_circles_lower_bound_K1

print("\n--- Final Calculation ---")
print("The difference between the braid index of K2 and the lower bound of the minimum number of Seifert circles of K1 is:")
print(f"{braid_index_K2} - {seifert_circles_lower_bound_K1} = {difference}")
print(f"<<<{difference}>>>")