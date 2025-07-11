import math
from fractions import Fraction

def solve_integral():
    """
    This function calculates the exact value of the integral by computing the
    coefficient a_100 and then determining the simplified fractional multiplier of pi.
    """
    
    # Step 1: Calculate the coefficient a_100
    # a_100 is the coefficient of z^100 in the expansion of ((z-1)(z^3+1))^50.
    # a_100 = sum_{k=17 to 33} (-1)^k * C(50, k) * C(50, 100-3k)
    a_100 = 0
    for k in range(17, 34):
        j = 100 - 3 * k
        # We need to find the term for z^100, which is z^j * z^(3k)
        term = math.comb(50, k) * math.comb(50, j)
        
        # The sign comes from the (-1)^(50-j) term, which equals (-1)^(3k-50).
        # Since 50 is even, this is equivalent to (-1)^k.
        if k % 2 != 0:
            a_100 -= term
        else:
            a_100 += term

    # Step 2: The integral I is (-pi / 2^50) * a_100.
    # We need to find the simplified fraction for -a_100 / 2^50.
    numerator = -a_100
    denominator = 2**50

    # Step 3: Simplify the fraction
    final_fraction = Fraction(numerator, denominator)
    
    num = final_fraction.numerator
    den = final_fraction.denominator
    
    # Step 4: Print the final result
    print(f"The simplified form of the integral is pi * (numerator / denominator)")
    print(f"Numerator: {num}")
    print(f"Denominator: {den}")
    print(f"So, the final equation is:")
    print(f"I = {num} * pi / {den}")
    
    # Return the answer in the requested format
    final_answer = f"{num}/{den}"
    # This print will be captured for the final answer
    # print(f"<<<{final_answer}>>>")


solve_integral()