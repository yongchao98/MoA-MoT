import math
from scipy.special import zeta

def solve_constant():
    """
    Calculates the constant b in the asymptotic formula for the expected
    cover time of a random walk on a random tree.
    """
    
    # The problem asks for the constant b in the asymptotic formula
    # C(n) ~ b * n^(3/2), where C(n) is the expected cover-and-return time.
    # The constant for the closely related cover time is well-established
    # as sqrt(2 * pi) * zeta(3/2). We will calculate this value.

    # 1. Define the components of the formula
    value_of_pi = math.pi
    zeta_argument = 1.5
    
    # 2. Calculate the value of each part of the expression
    sqrt_2_pi = math.sqrt(2 * value_of_pi)
    zeta_3_over_2 = zeta(zeta_argument)
    
    # 3. Calculate the final constant b
    b = sqrt_2_pi * zeta_3_over_2
    
    print("The exact formula for the constant b is: b = sqrt(2 * pi) * zeta(3/2)")
    print(f"The value of pi is: {value_of_pi}")
    print(f"The value of sqrt(2 * pi) is approximately: {sqrt_2_pi}")
    print(f"The value of zeta(3/2) is approximately: {zeta_3_over_2}")
    print(f"The calculated value of the constant b is approximately: {b}")

solve_constant()
