from fractions import Fraction

def solve():
    """
    Calculates the sum based on the analysis of the set S.
    The sum is S = (1/2) * [ 2*sum((1/4)^n) + sum((3/4)^n) + sum((-1/4)^n) ].
    """

    # For a geometric series sum_{n=1 to infinity} r^n, the sum is r / (1 - r).

    # Series 1: r = 1/4
    r1 = Fraction(1, 4)
    sum1 = r1 / (1 - r1)

    # Series 2: r = 3/4
    r2 = Fraction(3, 4)
    sum2 = r2 / (1 - r2)

    # Series 3: r = -1/4
    r3 = Fraction(-1, 4)
    sum3 = r3 / (1 - r3)

    # The total sum S = (1/2) * (2 * sum1 + sum2 + sum3)
    total_sum = Fraction(1, 2) * (2 * sum1 + sum2 + sum3)

    # Output the components of the final calculation
    print("The total sum is derived from three geometric series.")
    print(f"The sum of the series with r=1/4 is: {sum1}")
    print(f"The sum of the series with r=3/4 is: {sum2}")
    print(f"The sum of the series with r=-1/4 is: {sum3}")
    print("\nThe final value is calculated as:")
    
    # Print the equation with all the numbers
    term1 = 2 * sum1
    term2 = sum2
    term3 = sum3
    
    print(f"1/2 * (2 * {sum1} + {sum2} + {sum3}) = 1/2 * ({term1} + {term2} + ({term3})) = {total_sum}")

solve()