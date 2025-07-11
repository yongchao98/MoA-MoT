import math
from fractions import Fraction

solutions = []

def find_solutions_recursive(n, target_sum, min_x, current_nums):
    """
    Recursively finds sets of n distinct positive integers starting from min_x
    whose reciprocals sum to target_sum.
    """
    # Base case: if we only need to find one more number.
    if n == 1:
        # The last number x must satisfy 1/x = target_sum
        # This means x = 1/target_sum
        if target_sum.denominator == 1 and target_sum.numerator > 0:
            x = target_sum.numerator
            # The number must be distinct and larger than the previous one.
            if x > min_x:
                solutions.append(current_nums + [x])
        return

    # To ensure x_k is not too large or too small, we establish bounds.
    # Lower bound for the current number x: must be greater than the previous number.
    # Also, 1/x < target_sum / n (if all remaining numbers were equal to x).
    # So, x > n / target_sum
    lower_bound = max(min_x + 1, math.ceil(n / target_sum))
    
    # Upper bound for the current number x:
    # 1/x must be smaller than the target_sum itself.
    # So, x > 1 / target_sum.
    # And we also know target_sum = 1/x + ... < n/x, so x < n/target_sum is not an upper bound.
    # The actual upper bound comes from 1/x_k > target_sum / k.
    # It should be target_sum > n/x_{last} so x_{last} > n/target_sum
    # Let's use x_k < (k * x_{k-1}) / (k*...-1)
    # Simpler: 1/x < target_sum -> x > 1/target_sum
    # Let x_k be the current number to find. min_x is x_{k-1}.
    # We must have 1/x_k + ... + 1/x_n = target_sum. All x_i > min_x.
    # target_sum > n/x_n, so x_n > n/target_sum.
    # target_sum < n/x_k, so x_k < n/target_sum.
    upper_bound = math.floor(n / target_sum)

    for x in range(lower_bound, upper_bound + 1):
        new_target = target_sum - Fraction(1, x)
        if new_target >= 0:
            find_solutions_recursive(n - 1, new_target, x, current_nums + [x])

def solve_and_print():
    """
    Solves the T(4) problem and prints the results.
    """
    # Find all sets of 4 distinct integers {x1,x2,x3,x4} with 1 < x1 < x2 < x3 < x4
    # whose reciprocals sum to 1. Start with min_x = 1.
    find_solutions_recursive(4, Fraction(1), 1, [])
    
    print("Found the following sets of 4 distinct positive integers whose reciprocals sum to 1:\n")
    
    solution_sums = []
    for s in sorted(solutions, key=sum):
        s_str = ", ".join(map(str, s))
        reciprocal_str = " + ".join([f"1/{i}" for i in s])
        sum_val = sum(s)
        solution_sums.append(sum_val)
        
        print(f"Set: {{{s_str}}}")
        print(f"  - Equation: {reciprocal_str} = 1")
        print(f"  - Sum: {' + '.join(map(str, s))} = {sum_val}\n")

    print("-" * 40)
    
    # S(4) is the set of these sums.
    s4 = sorted(list(set(solution_sums)))
    s4_str = ", ".join(map(str, s4))
    print(f"The set of sums, S(4), is {{{s4_str}}}")
    
    # T(4) is the sum of the elements of S(4).
    t4 = sum(s4)
    t4_eq_str = " + ".join(map(str, s4))
    
    print("\nT(4) is the sum of these values.")
    print(f"Final Equation: T(4) = {t4_eq_str}")
    print(f"Result: T(4) = {t4}")

    print("\n<<<" + str(t4) + ">>>")

solve_and_print()