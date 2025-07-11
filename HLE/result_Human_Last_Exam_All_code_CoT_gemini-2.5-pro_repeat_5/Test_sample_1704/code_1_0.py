import math
from fractions import Fraction

def solve():
    """
    Finds T(4), the sum of all elements in S(4).
    S(n) is the set of sums of n distinct positive integers whose reciprocals sum to 1.
    """
    solutions = set()
    n = 4

    # The equation is 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1.
    # We assume 1 < x1 < x2 < x3 < x4.

    # Bound for x1: 1 = sum(1/xi) < n/x1 => x1 < n
    for x1 in range(2, n):
        rem1 = Fraction(1) - Fraction(1, x1)
        if rem1 <= 0:
            continue

        # Bound for x2: rem1 = sum(1/xi) < (n-1)/x2 => x2 < (n-1)/rem1
        x2_min = x1 + 1
        x2_max_limit = (n - 1) / rem1
        x2_max = math.ceil(x2_max_limit)
        
        for x2 in range(x2_min, x2_max):
            rem2 = rem1 - Fraction(1, x2)
            if rem2 <= 0:
                continue

            # Bound for x3: rem2 = sum(1/xi) < (n-2)/x3 => x3 < (n-2)/rem2
            # Also, x3 > 1/rem2 because rem2 = 1/x3 + 1/x4 > 1/x3
            x3_min = x2 + 1
            x3_lower_bound_from_rem = math.ceil(1 / rem2) if rem2 > 0 else x3_min
            x3_min = max(x3_min, x3_lower_bound_from_rem)

            x3_max_limit = (n - 2) / rem2 if rem2 > 0 else x3_min
            x3_max = math.ceil(x3_max_limit)

            for x3 in range(x3_min, x3_max):
                rem3 = rem2 - Fraction(1, x3)
                if rem3 <= 0:
                    continue
                
                # rem3 must be of the form 1/x4
                if rem3.numerator == 1:
                    x4 = rem3.denominator
                    if x4 > x3:
                        # Add the found set to our solutions
                        solutions.add(tuple(sorted((x1, x2, x3, x4))))

    # S(4) is the set of sums of the elements in each solution set.
    s4 = sorted([sum(sol) for sol in solutions])
    
    # T(4) is the sum of all elements in S(4).
    t4 = sum(s4)

    # Print the final equation as requested.
    equation_parts = [str(s) for s in s4]
    equation = " + ".join(equation_parts)
    print(f"{equation} = {t4}")

solve()