import math
from fractions import Fraction

def find_T4():
    """
    Finds T(4), the sum of all elements in S(4).
    S(n) is the set of all numbers that can be expressed as a sum of n distinct
    positive integers whose reciprocals sum to exactly 1.
    """
    
    solution_sets = []
    
    # Let the four distinct positive integers be x1, x2, x3, x4.
    # To find unique sets, we assume 1 < x1 < x2 < x3 < x4.
    # From 1 = 1/x1 + 1/x2 + 1/x3 + 1/x4 < 4/x1, we get x1 < 4.
    # x1 cannot be 1, so x1 can be 2 or 3.
    
    # Loop for x1
    for x1 in range(2, 4):
        
        # Loop for x2
        # From 1 - 1/x1 = 1/x2 + 1/x3 + 1/x4 < 3/x2, we can derive a limit for x2.
        # x2 < 3 / (1 - 1/x1)
        limit_x2 = math.ceil(3 / (Fraction(1) - Fraction(1, x1)))
        for x2 in range(x1 + 1, limit_x2):
            
            # Loop for x3
            # From 1 - 1/x1 - 1/x2 = 1/x3 + 1/x4 < 2/x3, we derive a limit for x3.
            # x3 < 2 / (1 - 1/x1 - 1/x2)
            rem_x3 = Fraction(1) - Fraction(1, x1) - Fraction(1, x2)
            if rem_x3 <= 0:
                continue
            
            # We also need 1/x3 < rem_x3, which means x3 > 1/rem_x3
            start_x3 = math.floor(1 / rem_x3) + 1
            limit_x3 = math.ceil(2 / rem_x3)
            
            for x3 in range(max(x2 + 1, start_x3), limit_x3):
                
                # Calculate the required 1/x4
                rem_x4 = Fraction(1) - Fraction(1, x1) - Fraction(1, x2) - Fraction(1, x3)
                
                # Check if it's a valid unit fraction for a new distinct integer
                if rem_x4 > 0 and rem_x4.numerator == 1:
                    x4 = rem_x4.denominator
                    if x4 > x3:
                        solution_sets.append((x1, x2, x3, x4))

    print("Finding T(4), the sum of all elements in S(4).")
    print("S(4) is the set of sums of 4 distinct positive integers whose reciprocals sum to 1.")
    print("-" * 30)

    total_sum_T4 = 0
    sum_elements = []

    # Sort for a consistent output order
    solution_sets.sort()
    
    for s in solution_sets:
        x1, x2, x3, x4 = s
        current_sum = sum(s)
        sum_elements.append(current_sum)
        total_sum_T4 += current_sum
        
        print(f"Found solution set: {{{x1}, {x2}, {x3}, {x4}}}")
        print(f"Equation: 1/{x1} + 1/{x2} + 1/{x3} + 1/{x4} = 1")
        print(f"Sum: {x1} + {x2} + {x3} + {x4} = {current_sum}")
        print("-" * 30)
        
    print(f"The set of sums, S(4), is: {sorted(sum_elements)}")
    print(f"T(4) is the sum of these elements:")
    print(f"{' + '.join(map(str, sorted(sum_elements)))} = {total_sum_T4}")
    
    # Return the final answer in the specified format
    return total_sum_T4

if __name__ == '__main__':
    final_answer = find_T4()
    print(f"\n<<< {final_answer} >>>")