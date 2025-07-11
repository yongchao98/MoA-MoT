import math
from fractions import Fraction

def compute_sum_coefficient():
    """
    Computes the rational coefficient of the sum L = zeta(6)*zeta(8)/zeta(24).
    The sum is of the form C * pi^(-10).
    This function computes and prints C.
    """
    
    # Values of Bernoulli numbers
    # |B_6| = 1/42
    # |B_8| = 1/30
    # |B_24| = 236364091/2730
    
    # We will compute the coefficient C = (zeta(6)/pi^6) * (zeta(8)/pi^8) / (zeta(24)/pi^24)
    # C = (1/945) * (1/9450) / ( |B_24| * 2^23 / 24! )
    
    B24_num = 236364091
    B24_den = 2730

    # Denominator of zeta(24)/pi^24: |B_24|*2^23 / 24!
    zeta24_coeff_num = Fraction(B24_num, B24_den) * Fraction(2**23)
    zeta24_coeff_den = Fraction(math.factorial(24))
    
    zeta24_coeff = zeta24_coeff_num / zeta24_coeff_den

    # Coefficient for zeta(6)/pi^6 is 1/945
    zeta6_coeff = Fraction(1, 945)
    
    # Coefficient for zeta(8)/pi^8 is 1/9450
    zeta8_coeff = Fraction(1, 9450)

    # Total rational coefficient C
    C = (zeta6_coeff * zeta8_coeff) / zeta24_coeff
    
    print(f"The sum is given by the expression zeta(6)*zeta(8)/zeta(24).")
    print(f"zeta(6) = (pi^6)/945")
    print(f"zeta(8) = (pi^8)/9450")
    print(f"zeta(24) = ({B24_num}/{B24_den}) * (2^23 * pi^24) / (24!)")
    print(f"The sum evaluates to C * pi^(-10), where C is:")
    print(f"C = {C.numerator}/{C.denominator}")
    
    final_expression = f"{C.numerator}/{C.denominator} * pi^-10"
    print("\nFinal Answer:")
    print(final_expression)

compute_sum_coefficient()