import sympy
from pyknotid.representations import Braid
import pyknotid.catalogue as cat

def get_conway_z2_coefficient(alexander_poly, t_var):
    """
    Calculates the z^2 coefficient of the Conway polynomial from the Alexander polynomial.
    
    The Alexander polynomial P(t) from pyknotid is a polynomial in t.
    It corresponds to a symmetric Laurent polynomial Delta(t).
    This function converts P(t) to the normalized Delta(t) and calculates the
    z^2 coefficient of the corresponding Conway polynomial.
    
    The formula for the z^2 coefficient is sum(j^2 * a_j) where a_j are the
    coefficients of the normalized Laurent polynomial Delta(t) = sum(a_j * t^j).
    """
    if not isinstance(alexander_poly, sympy.Expr) or alexander_poly.is_number:
        return 0

    p = sympy.Poly(alexander_poly, t_var)
    
    # The Alexander polynomial of a knot, evaluated at t=1, is +/- 1.
    # We normalize it to be 1, as required for the Conway polynomial correspondence.
    val_at_1 = p.eval(1)
    if val_at_1 == 0:
        # This would happen for a link with multiple components.
        # The braid closure is a knot, so this case shouldn't be reached.
        raise ValueError("Alexander polynomial is 0 at t=1.")

    normalized_poly = p / val_at_1
    
    # The degree of the polynomial from pyknotid. It's symmetric.
    degree = normalized_poly.degree()
    if degree % 2 != 0:
        # This implies the center of symmetry of the Laurent polynomial is not t^0,
        # but pyknotid's standard output for knots is a symmetric polynomial with even degree.
        raise ValueError("Polynomial has odd degree, cannot form a standard symmetric Laurent polynomial.")
        
    d_half = degree // 2
    
    # Calculate the z^2 coefficient using the formula sum(j^2 * a_j).
    # a_j is the coefficient of t^j in the Laurent polynomial Delta(t).
    # Delta(t) = t^(-d_half) * normalized_poly(t).
    # So, a_j = coefficient of t^(j+d_half) in normalized_poly.
    
    coeff_z2 = 0
    for j in range(1, d_half + 1):
        k = j + d_half
        a_j = normalized_poly.coeff_monomial(t_var**k)
        coeff_z2 += j*j * a_j
        
    return coeff_z2

# Define the symbolic variable for our polynomials
t = sympy.Symbol('t')

# --- Part 1: Analysis of the braid beta ---

# Define the braid beta from the problem statement
beta_word = [-4, -4, -3, 4, -3, 2, -1, -3, -2, -2, -2, -1]
# This is a braid on 5 strands
beta_braid = Braid(beta_word, num_strands=5)

# Calculate the Alexander polynomial for the closure of beta
poly_beta = beta_braid.alexander_polynomial(variable=t)

# Calculate the z^2 coefficient of the Conway polynomial for the braid closure
coeff_beta = get_conway_z2_coefficient(poly_beta, t)


# --- Part 2: Analysis of the knot 10_4 ---

# Get the knot 10_4 from the pyknotid catalogue
k10_4, = cat.get_knot('10_4')

# Get the Alexander polynomial for 10_4
poly_10_4 = k10_4.alexander_polynomial(variable=t)

# Calculate the z^2 coefficient of the Conway polynomial for 10_4
coeff_10_4 = get_conway_z2_coefficient(poly_10_4, t)


# --- Part 3: Calculate the difference and print the results ---

difference = coeff_beta - coeff_10_4

print(f"The Alexander polynomial for the closure of beta is: {poly_beta}")
print(f"The Alexander polynomial for the knot 10_4 is: {poly_10_4}")
print("")
print("The z^2 coefficient of the Conway polynomial for the closure of beta is:", coeff_beta)
print("The z^2 coefficient of the Conway polynomial for the knot 10_4 is:", coeff_10_4)
print("")
print(f"The difference is: {coeff_beta} - ({coeff_10_4}) = {difference}")
