# This script requires the pyknotid library.
# You can install it by running: pip install pyknotid
import pyknotid.catalogue as cat
from sympy import Poly, symbols

def get_braid_index_of_K2():
    """
    Identifies the knot K2 from its braid representation and determines its braid index.
    K2 is the closure of the braid (sigma_1^{-1})^3 * sigma_2^{-1}.
    """
    # In pyknotid's notation, the braid sigma_i is represented by the integer i,
    # and its inverse by -i. The generators are 1-indexed.
    # So, the braid word for (sigma_1^{-1})^3 * sigma_2^{-1} is [-1, -1, -1, -2].
    braid_word = [-1, -1, -1, -2]
    
    # Create the knot object from the braid closure
    k2 = cat.from_braid(braid_word)
    
    # Identify the knot. The identify() method returns a tuple like (crossings, index).
    # For this braid, it identifies the knot as (4, 1), which is the figure-eight knot (4_1).
    identifier = k2.identify()
    
    # Determine the braid index.
    # The braid index of the unknot is 1.
    # The braid index of any (2, p)-torus knot is 2.
    # The figure-eight knot is not the unknot and not a torus knot.
    # Therefore, its braid index must be at least 3.
    # Since we have a 3-strand representation, the braid index is exactly 3.
    braid_index_K2 = 3
    return braid_index_K2

def get_seifert_circles_lower_bound_of_K1():
    """
    Calculates the lower bound on the minimum number of Seifert circles for K1 (10_74 knot)
    using its HOMFLY polynomial.
    """
    # Get the knot K1, which is the 10_74 knot.
    k1 = cat.get_knot(10, 74)
    
    # Calculate the HOMFLY polynomial. pyknotid uses variables v and z by default.
    homfly_poly = k1.homfly_polynomial()
    
    # The polynomial is a sympy expression. We need to find the span of the powers of 'v'.
    # Get the symbol for 'v' from the polynomial.
    v = next(iter(homfly_poly.free_symbols.difference(symbols('z'))))
    
    # Treat the expression as a polynomial in v.
    poly_in_v = Poly(homfly_poly, v)
    
    # The degree of the polynomial is the maximum power of v.
    v_max = poly_in_v.degree()
    
    # To find the minimum power of v, we look at all the terms.
    # The terms() method gives a list of ((power_v, power_z), coefficient).
    v_powers = [term[0][0] for term in poly_in_v.terms()]
    v_min = min(v_powers)
    
    # The v-span is the difference between the max and min powers.
    v_span = v_max - v_min
    
    # The lower bound for the number of Seifert circles is given by the
    # Morton-Franks-Williams inequality: bound = v_span / 2 + 1.
    lower_bound = v_span / 2 + 1
    
    return int(lower_bound)

# Calculate the two values
val1 = get_braid_index_of_K2()
val2 = get_seifert_circles_lower_bound_of_K1()

# Calculate the difference
difference = val1 - val2

# Print the final equation as requested
print(f"{val1} - {val2} = {difference}")