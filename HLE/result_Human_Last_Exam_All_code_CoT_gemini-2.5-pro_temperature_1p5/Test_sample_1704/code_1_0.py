import math
from fractions import Fraction

def find_t4():
    """
    Finds all sets of 4 distinct positive integers whose reciprocals sum to 1,
    calculates the sum of each set (these form the set S(4)),
    and then calculates the sum of the elements of S(4) (which is T(4)).
    """
    n = 4
    solutions = set()

    # We are solving 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1
    # Assume 1 < x1 < x2 < x3 < x4.
    # From analysis, x1 must be 2 or 3.
    for x1 in range(2, 4):
        rem1 = Fraction(1) - Fraction(1, x1)
        
        # Determine upper bound for x2: x2 < (n-1)/rem1
        limit2 = (n - 1) / rem1
        x2_upper = math.ceil(float(limit2))
        for x2 in range(x1 + 1, x2_upper):
            rem2 = rem1 - Fraction(1, x2)
            if rem2 <= 0:
                continue

            # Determine upper bound for x3: x3 < (n-2)/rem2
            limit3 = (n - 2) / rem2
            x3_upper = math.ceil(float(limit3))
            for x3 in range(x2 + 1, x3_upper):
                rem3 = rem2 - Fraction(1, x3)
                
                # rem3 must be 1/x4, so its numerator must be 1
                if rem3 > 0 and rem3.numerator == 1:
                    x4 = rem3.denominator
                    # Ensure the integers are distinct
                    if x4 > x3:
                        solutions.add(tuple(sorted((x1, x2, x3, x4))))

    # S(4) is the set of sums of the integers in each solution.
    s4_sums = {sum(s) for s in solutions}

    # T(4) is the sum of the elements in S(4).
    t4_total = sum(s4_sums)

    # Sort the sums for a clean presentation.
    sorted_s4 = sorted(list(s4_sums))

    # Build and print the final equation as requested.
    equation_str = " + ".join(map(str, sorted_s4))
    print(f"The elements of S(4) are {sorted_s4}.")
    print(f"T(4) is the sum of these elements:")
    print(f"{equation_str} = {t4_total}")

find_t4()