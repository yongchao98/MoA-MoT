import math

def solve():
    """
    Calculates the integer part of 10^4 * L, where L is the limit described in the problem.
    """
    # For the family of solutions with k=13
    k1 = 13
    # The growth rate is given by the larger root of the characteristic equation x^2 - kx + 1 = 0
    lambda1 = (k1 + math.sqrt(k1**2 - 4)) / 2
    # The contribution to the limit from this family
    limit_term1 = 2 / math.log(lambda1)

    # For the family of solutions with k=5
    k2 = 5
    # The growth rate for this family
    lambda2 = (k2 + math.sqrt(k2**2 - 4)) / 2
    # The contribution to the limit from this family
    limit_term2 = 2 / math.log(lambda2)

    # The total limit L is the sum of the contributions from the two families.
    L = limit_term1 + limit_term2

    # We need to find the integer part of 10^4 * L.
    final_value = math.floor(10000 * L)

    # Outputting the numbers in the final equation as requested.
    print("The final limit L is determined by the equation:")
    print(f"L = 2 / ln(({k1} + sqrt({k1**2 - 4})) / 2) + 2 / ln(({k2} + sqrt({k2**2 - 4})) / 2)")
    print(f"\nThis simplifies to:")
    print(f"L = 2 / ln(({k1} + {math.sqrt(k1**2 - 4)}) / 2) + 2 / ln(({k2} + {math.sqrt(k2**2 - 4)}) / 2)")
    print(f"L = 2 / {math.log(lambda1)} + 2 / {math.log(lambda2)}")
    print(f"L = {limit_term1} + {limit_term2}")
    print(f"L = {L}")
    print(f"\nThe value of 10^4 * L is {10000 * L}")
    print(f"The integer part of 10^4 * L is: {final_value}")

solve()