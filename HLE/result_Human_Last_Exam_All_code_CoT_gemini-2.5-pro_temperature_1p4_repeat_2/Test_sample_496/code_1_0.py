import math

def calculate_total_rank():
    """
    Calculates the total rank of the SO(4)-equivariant cohomology ring of SO(4)\X
    in degrees up to 100.

    The Poincare series of the cohomology ring is given by P(t) = (t^3 + t^6) / (1 - t^4)^2.
    This can be expanded as:
    P(t) = (t^3 + t^6) * sum_{k=0 to infinity} (k+1) * t^(4k)
         = sum_{k=0 to infinity} (k+1) * t^(4k + 3) + sum_{k=0 to infinity} (k+1) * t^(4k + 6)

    We sum the coefficients (ranks) for all terms with degree at most 100.
    """

    degree_limit = 100
    total_rank = 0

    # Contribution from the first series: sum (k+1)t^(4k+3)
    # We need 4k + 3 <= 100, which implies 4k <= 97, so k <= 24.25.
    # k ranges from 0 to 24.
    max_k1 = math.floor((degree_limit - 3) / 4)
    sum1 = 0
    # The coefficients are k+1, for k = 0, ..., max_k1
    # This is the sum of integers from 1 to max_k1 + 1
    if max_k1 >= 0:
        sum1 = (max_k1 + 1) * (max_k1 + 2) // 2
    
    # Contribution from the second series: sum (k+1)t^(4k+6)
    # We need 4k + 6 <= 100, which implies 4k <= 94, so k <= 23.5.
    # k ranges from 0 to 23.
    max_k2 = math.floor((degree_limit - 6) / 4)
    sum2 = 0
    # The coefficients are k+1, for k = 0, ..., max_k2
    # This is the sum of integers from 1 to max_k2 + 1
    if max_k2 >= 0:
        sum2 = (max_k2 + 1) * (max_k2 + 2) // 2

    total_rank = sum1 + sum2
    
    print(f"The sum of ranks from the first part of the formula is: {sum1}")
    print(f"The sum of ranks from the second part of the formula is: {sum2}")
    print(f"The total rank is the sum of these two numbers.")
    print(f"{sum1} + {sum2} = {total_rank}")

calculate_total_rank()