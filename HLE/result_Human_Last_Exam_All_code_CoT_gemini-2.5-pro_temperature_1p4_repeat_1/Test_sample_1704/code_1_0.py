import math
from fractions import Fraction

def find_t4():
    """
    Finds the value of T(4), where T(n) is the sum of all elements in S(n),
    and S(n) is the set of all numbers that can be expressed as a sum of n
    distinct positive integers whose reciprocals sum to exactly 1.
    """
    n = 4
    solutions = []

    # We are looking for distinct positive integers x1, x2, x3, x4 such that
    # 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1.
    # To avoid duplicates, we enforce an order: 1 < x1 < x2 < x3 < x4.

    # Bound for x1: 1 = 1/x1 + 1/x2 + 1/x3 + 1/x4 < 4/x1 => x1 < 4.
    # So, x1 can be 2 or 3.
    for x1 in range(2, n):
        
        # Bound for x2: 1 - 1/x1 = 1/x2 + 1/x3 + 1/x4 < 3/x2 => x2 < 3 / (1 - 1/x1).
        x2_upper_bound = 3 / (Fraction(1) - Fraction(1, x1))
        for x2 in range(x1 + 1, math.ceil(x2_upper_bound)):
            
            # Bound for x3: 1 - 1/x1 - 1/x2 = 1/x3 + 1/x4 < 2/x3 => x3 < 2 / (1 - 1/x1 - 1/x2).
            rem_for_x3_x4 = Fraction(1) - Fraction(1, x1) - Fraction(1, x2)
            if rem_for_x3_x4 <= 0:
                continue
            x3_upper_bound = 2 / rem_for_x3_x4
            for x3 in range(x2 + 1, math.ceil(x3_upper_bound)):
                
                # Now, we calculate the required 1/x4.
                rem_for_x4 = Fraction(1) - Fraction(1, x1) - Fraction(1, x2) - Fraction(1, x3)

                # For x4 to be a distinct positive integer, rem_for_x4 must be a positive
                # unit fraction (numerator is 1), and its denominator must be > x3.
                if rem_for_x4 > 0 and rem_for_x4.numerator == 1:
                    x4 = rem_for_x4.denominator
                    if x4 > x3:
                        solution = (x1, x2, x3, x4)
                        solutions.append(solution)
    
    print("Found all sets of 4 distinct positive integers whose reciprocals sum to 1:\n")
    
    s4 = set()
    for sol in sorted(solutions):
        s = sum(sol)
        s4.add(s)
        print(f"Set: {sol}")
        print(f"  - Reciprocal sum: 1/{sol[0]} + 1/{sol[1]} + 1/{sol[2]} + 1/{sol[3]} = 1")
        print(f"  - Sum of elements (element of S(4)): {sol[0]} + {sol[1]} + {sol[2]} + {sol[3]} = {s}\n")

    sorted_s4 = sorted(list(s4))
    t4 = sum(sorted_s4)

    print("-" * 30)
    print(f"The set S(4) is the set of these sums: {sorted_s4}")
    
    sum_equation = " + ".join(map(str, sorted_s4))
    print(f"T(4) is the sum of all elements in S(4):")
    print(f"T(4) = {sum_equation} = {t4}")

if __name__ == "__main__":
    find_t4()
    print("\n<<<208>>>")
