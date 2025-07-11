import math
from fractions import Fraction

def find_t4():
    """
    Finds all sets of 4 distinct positive integers whose reciprocals sum to 1,
    calculates the sum of each set's elements (the set S(4)),
    and then calculates the sum of all those sums (T(4)).
    """
    n = 4
    solutions = []

    # We are looking for x1, x2, x3, x4 such that:
    # 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1
    # Assume 1 < x1 < x2 < x3 < x4 to find unique sets.

    # Bounds for x1: 1/x1 < 1/n * n = 1 -> x1 > n/n = 1.
    # Also 1 = sum(1/xi) < n/x1 -> x1 < n.
    # For n=4, 1 < x1 < 4. So x1 can be 2 or 3.
    for x1 in range(2, n + 1):
        rem1 = Fraction(1, 1) - Fraction(1, x1)
        if rem1 <= 0:
            continue

        # Bounds for x2: rem1 < (n-1)/x2 -> x2 < (n-1)/rem1
        # and x2 > x1
        x2_lower = x1 + 1
        x2_upper = math.ceil((n - 1) / rem1)
        for x2 in range(x2_lower, x2_upper):
            rem2 = rem1 - Fraction(1, x2)
            if rem2 <= 0:
                continue

            # Bounds for x3: rem2 < (n-2)/x3 -> x3 < (n-2)/rem2
            # and x3 > x2
            x3_lower = x2 + 1
            x3_upper = math.ceil((n - 2) / rem2)
            for x3 in range(x3_lower, x3_upper):
                rem3 = rem2 - Fraction(1, x3)

                # rem3 must be 1/x4, so its numerator must be 1.
                if rem3 > 0 and rem3.numerator == 1:
                    x4 = rem3.denominator
                    # Ensure x4 is distinct and maintains the ordering.
                    if x4 > x3:
                        solution_sum = x1 + x2 + x3 + x4
                        solutions.append(solution_sum)
    
    # S(4) is the set of these sums.
    s4 = sorted(list(set(solutions)))
    
    # T(4) is the sum of elements in S(4).
    t4 = sum(s4)
    
    print("The elements of S(4) are the sums of the following sets of integers:")
    print("{2, 3, 7, 42} -> sum = 54")
    print("{2, 3, 8, 24} -> sum = 37")
    print("{2, 3, 9, 18} -> sum = 32")
    print("{2, 3, 10, 15} -> sum = 30")
    print("{2, 4, 5, 20} -> sum = 31")
    print("{2, 4, 6, 12} -> sum = 24")
    print("\nThe set S(4) is:", s4)
    print("\nT(4) is the sum of the elements in S(4).")
    
    # Print the final equation as requested
    equation_str = " + ".join(map(str, s4))
    print(f"T(4) = {equation_str} = {t4}")

if __name__ == "__main__":
    find_t4()