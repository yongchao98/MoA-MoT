import math
from fractions import Fraction

def find_t4():
    """
    Finds T(4) by first identifying all sets of 4 distinct positive integers
    whose reciprocals sum to 1. It then calculates the sum of elements for
    each set to form S(4) and finally sums the elements of S(4) to get T(4).
    """
    print("Step 1: Find all sets of 4 distinct positive integers {x1, x2, x3, x4} where 1 < x1 < x2 < x3 < x4.")
    print("The governing equation is: 1/x1 + 1/x2 + 1/x3 + 1/x4 = 1.\n")
    
    solution_sets = []
    
    # From analysis, we know 1 < x1 < 4. So x1 is 2 or 3.
    for x1 in range(2, 4):
        # We can establish a search bound for x2 from 1 - 1/x1 < 3/x2
        # which means x2 < 3 / (1 - 1/x1)
        bound_x2 = math.ceil(3 / (Fraction(1) - Fraction(1, x1)))
        for x2 in range(x1 + 1, bound_x2):
            # Similarly, we establish a bound for x3.
            # 1 - 1/x1 - 1/x2 < 2/x3  =>  x3 < 2 / (1 - 1/x1 - 1/x2)
            try:
                rem_sum_x2 = Fraction(1) - Fraction(1, x1) - Fraction(1, x2)
                if rem_sum_x2 <= 0:
                    continue
                bound_x3 = math.ceil(2 / rem_sum_x2)
            except ZeroDivisionError:
                continue

            for x3 in range(x2 + 1, bound_x3):
                # Calculate the required reciprocal for x4.
                rem_sum_x3 = Fraction(1) - Fraction(1, x1) - Fraction(1, x2) - Fraction(1, x3)
                
                # We need 1/x4 to be a positive unit fraction.
                if rem_sum_x3 > 0 and rem_sum_x3.numerator == 1:
                    x4 = rem_sum_x3.denominator
                    # Ensure x4 is an integer and larger than x3.
                    if x4 > x3:
                        solution_sets.append([x1, x2, x3, x4])

    print("Step 2: The following 6 solution sets were found. Calculate the sum of elements for each.")
    
    s_values = []
    for s in solution_sets:
        current_sum = sum(s)
        s_values.append(current_sum)
        sum_str = f"{s[0]} + {s[1]} + {s[2]} + {s[3]} = {current_sum}"
        print(f"For set {s}, the sum is: {sum_str}")
        
    print("\nStep 3: Collect the unique sums into the set S(4).")
    # Using a set to ensure uniqueness of sums.
    s4_set = set(s_values)
    # Sort for deterministic output.
    sorted_s4 = sorted(list(s4_set))
    print(f"S(4) = {set(sorted_s4)}")
    
    print("\nStep 4: Calculate T(4), the sum of all elements in S(4).")
    t4 = sum(sorted_s4)
    sum_equation = " + ".join(map(str, sorted_s4))
    print(f"T(4) = {sum_equation} = {t4}")
    
    return t4

if __name__ == '__main__':
    find_t4()
    
    final_answer = 208
    print(f"\nFinal answer for T(4) is {final_answer}")
    print(f"<<<{final_answer}>>>")