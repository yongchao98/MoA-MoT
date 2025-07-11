def solve_problem():
    """
    This function calculates the final value based on the derived formula.
    The problem simplifies to calculating the sum of the first n integers, where n=4048,
    and adding 10^15 to the result.
    """
    n = 4048
    
    # Calculate the sum of alpha_i^2, which is sum of i from 1 to n
    # This is given by the formula n * (n + 1) / 2
    sum_of_alphas_squared = n * (n + 1) // 2
    
    # The second term in the final expression
    large_number = 10**15
    
    # Calculate the final result
    result = sum_of_alphas_squared + large_number
    
    # Print the equation with each number as requested
    print(f"{sum_of_alphas_squared} + {large_number} = {result}")

solve_problem()