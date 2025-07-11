import math

def solve():
    """
    This function calculates the final value based on the derived formula.
    """
    n = 4048
    # The value to be calculated simplifies to sum(alpha_i^2) + 10^15.
    # Given alpha_i = sqrt(i), alpha_i^2 = i.
    # We need to calculate sum_{i=1 to n} i.
    
    # Calculate the sum of the first n integers
    # The result must be an integer as either n or n+1 is even.
    sum_alpha_squared = n * (n + 1) // 2
    
    # The constant term
    constant_term = 10**15
    
    # Calculate the final result
    final_result = sum_alpha_squared + constant_term
    
    # Print the equation with all the numbers, as requested.
    print(f"{sum_alpha_squared} + {constant_term} = {final_result}")

solve()