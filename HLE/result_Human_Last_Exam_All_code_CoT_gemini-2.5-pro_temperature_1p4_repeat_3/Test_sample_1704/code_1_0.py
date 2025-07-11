import math
from fractions import Fraction

def find_t4():
    """
    This function finds all sets of 4 distinct positive integers {x1, x2, x3, x4}
    such that their reciprocals sum to 1. It then calculates S(4), the set of
    the sums of these integers, and T(4), the sum of the elements in S(4).
    """
    n = 4
    solutions = []
    
    # We are solving 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1
    # Assume x1 < x2 < x3 < x4 to avoid duplicates and simplify search.
    
    # Bound for x1: 1 < x1 < n
    # The loop for x1 will go from 2 to 3.
    for x1 in range(2, n):
        # Bound for x2: x1 < x2 < (n-1)*x1 / (x1 - 1)
        x2_upper_bound = math.ceil((n - 1) * x1 / (x1 - 1.0))
        for x2 in range(x1 + 1, x2_upper_bound):
            # Remaining value for the rest of the terms
            rem1 = Fraction(1) - Fraction(1, x1) - Fraction(1, x2)
            if rem1 <= 0:
                continue

            # Bound for x3: x2 < x3 < (n-2) / rem1
            x3_upper_bound = math.ceil((n - 2) / rem1)
            # Also x3 > 1/rem1
            x3_lower_bound = math.floor(1 / rem1) + 1
            
            for x3 in range(max(x2 + 1, x3_lower_bound), x3_upper_bound):
                rem2 = rem1 - Fraction(1, x3)
                if rem2 <= 0:
                    continue
                
                # We need rem2 = 1/x4. This means rem2.numerator must be 1.
                if rem2.numerator == 1:
                    x4 = rem2.denominator
                    # Ensure x4 is distinct and maintains the order.
                    if x4 > x3:
                        solutions.append((x1, x2, x3, x4))

    # S(4) is the set of sums of the solution sets.
    s4_sums = set()
    for sol in solutions:
        s4_sums.add(sum(sol))
    
    # Sort for consistent and readable output
    sorted_sums = sorted(list(s4_sums))
    
    print("This script finds T(4) by first identifying all sets of 4 distinct positive integers whose reciprocals sum to 1.")
    print("The equation is: 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1\n")
    print("The found sets {x1, x2, x3, x4} and their corresponding sums are:")
    
    # Sort solutions by their sum for ordered printing
    sorted_solutions = sorted(solutions, key=lambda x: sum(x))
    
    for sol in sorted_solutions:
        s = sum(sol)
        print(f"Set: {sol}, Sum: {sol[0]} + {sol[1]} + {sol[2]} + {sol[3]} = {s}")

    print("\nThe set of unique sums, S(4), is:")
    print(s4_sums)

    t4 = sum(s4_sums)
    print("\nT(4) is the sum of these unique values:")
    sum_str = " + ".join(map(str, sorted_sums))
    print(f"T(4) = {sum_str} = {t4}")

if __name__ == '__main__':
    find_t4()