import math
from fractions import Fraction

def solve():
    """
    This function calculates the rational coefficient for the sum Σ = ζ(6)ζ(8)/ζ(24).
    The sum is expressed as (coefficient) * π^(-10).
    """
    
    # The sum Σ can be expressed as ζ(6)ζ(8)/ζ(24).
    # We use the formula ζ(2k) = |B_2k| * (2π)^(2k) / (2 * (2k)!).
    # This leads to Σ = R * π^(-10), where R is a rational number.
    # R = (c_6 * c_8) / c_24, where c_{2k} = |B_{2k}| * 2^(2k-1) / (2k)!
    # R = (|B_6|*|B_8|/|B_{24}|) * (2^(5+7-23)) * (24!/(6!*8!))
    # R = (|B_6|*|B_8|/|B_{24}|) * 2^(-11) * (24!/(6!*8!))

    # Known absolute values of Bernoulli numbers as Fractions for precision
    # B_6 = 1/42
    # B_8 = -1/30
    # B_24 = -236364091/2730
    B6_abs = Fraction(1, 42)
    B8_abs = Fraction(1, 30)
    B24_abs = Fraction(236364091, 2730)

    # Calculate the ratio of Bernoulli numbers
    B_ratio = (B6_abs * B8_abs) / B24_abs

    # Calculate the power of 2 term
    pow2_term = Fraction(1, 2**11)

    # Calculate the factorial term
    fact_term = Fraction(math.factorial(24), math.factorial(6) * math.factorial(8))

    # Combine all parts to get the final rational coefficient
    rational_coefficient = B_ratio * pow2_term * fact_term

    num = rational_coefficient.numerator
    den = rational_coefficient.denominator
    
    print("The sum is given by the expression: zeta(6) * zeta(8) / zeta(24)")
    print(f"This evaluates to the final answer:")
    print(f"({num} / {den}) * pi^(-10)")

solve()
<<<22610499275/236364091 * pi**-10>>>