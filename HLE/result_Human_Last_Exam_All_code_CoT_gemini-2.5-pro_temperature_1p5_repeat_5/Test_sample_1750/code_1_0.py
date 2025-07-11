import math

def calculate_integral_fraction():
    """
    This function calculates the rational part of the definite integral.
    The integral simplifies to (pi * q), and this function finds q.
    The value q is -C_100 / 2^50, where C_100 is a sum of binomial coefficients.
    """
    
    # C_100 is the coefficient of z^100 in ((z^3+1)(z-1))^50
    # C_100 = sum_{j=17 to 33} (-1)^j * comb(50, j) * comb(50, 100-3j)
    C100 = 0
    for j in range(17, 34):
        # The combinatorial identity is Sum[(-1)^j * C(n,j) * C(n,k-cj)]
        # I have k = 100, n=50, c=3.
        # (-1)^(3j-50) is the exponent on -1. Since 50 is even, this is (-1)^(3j) = (-1)^j.
        term = math.comb(50, j) * math.comb(50, 100 - 3 * j)
        if j % 2 != 0:
            C100 -= term
        else:
            C100 += term

    numerator = -C100
    denominator = 2**50
    
    # Simplify the fraction
    common_divisor = math.gcd(numerator, denominator)
    
    num_simplified = numerator // common_divisor
    den_simplified = denominator // common_divisor
    
    print("The integral is equal to the following fraction multiplied by pi:")
    print(f"Numerator: {num_simplified}")
    print(f"Denominator: {den_simplified}")
    print("\nThe final equation is:")
    print(f"Integral = {num_simplified}/{den_simplified} * pi")

calculate_integral_fraction()
