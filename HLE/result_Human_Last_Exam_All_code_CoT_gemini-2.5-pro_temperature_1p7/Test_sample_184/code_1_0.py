import math
from fractions import Fraction

def bernoulli(n):
    """
    Returns the n-th Bernoulli number.
    This function is only implemented for the required values.
    B_6 = 1/42
    B_8 = -1/30
    B_24 = 236364091/2730
    """
    if n == 6:
        return Fraction(1, 42)
    elif n == 8:
        # The formula for zeta(2k) uses the absolute value |B_2k|
        return Fraction(1, 30)
    elif n == 24:
        # The formula for zeta(2k) uses the absolute value |B_2k|
        return Fraction(236364091, 2730)
    else:
        raise NotImplementedError("Bernoulli number not implemented for n=" + str(n))

def zeta_rational_part(n):
    """
    Calculates zeta(n) / pi^n, which is a rational number for even n.
    zeta(2k) / pi^(2k) = 2^(2k-1) * |B_2k| / (2k)!
    """
    if n % 2 != 0:
        raise ValueError("Can only compute for even n")
    k = n // 2
    abs_B_2k = bernoulli(n)
    numerator = 2**(n - 1) * abs_B_2k.numerator
    denominator = math.factorial(n) * abs_B_2k.denominator
    return Fraction(numerator, denominator)

# The sum is zeta(6)*zeta(8)/zeta(24)
# We calculate the rational coefficient C for C * pi^(-10)
zeta_star_6 = zeta_rational_part(6)
zeta_star_8 = zeta_rational_part(8)
zeta_star_24 = zeta_rational_part(24)

# C = (zeta(6)/pi^6) * (zeta(8)/pi^8) / (zeta(24)/pi^24)
rational_coeff = zeta_star_6 * zeta_star_8 / zeta_star_24

num = rational_coeff.numerator
den = rational_coeff.denominator

# The final answer is C * pi^(-10)
print("The sum is given by the expression: zeta(6) * zeta(8) / zeta(24)")
print(f"zeta(6) = pi^6 / {zeta_star_6.denominator}")
print(f"zeta(8) = pi^8 / {zeta_star_8.denominator}")
# The expression for zeta(24) is complicated, so we just show the result.
print(f"The final sum evaluates to ({num} / {den}) * pi^-10")
