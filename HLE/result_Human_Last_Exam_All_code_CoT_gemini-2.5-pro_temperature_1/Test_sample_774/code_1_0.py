import numpy as np

def solve_problem():
    """
    This function calculates the final result based on the step-by-step analysis.
    """
    
    # Part 1: The Summation
    # The term l(n,p) is the injectivity radius of the Stiefel manifold M(n,p), which is pi.
    # The sum is over a 10x10 grid, resulting in 100 * pi.
    sum_coefficient = 100
    pi_value = np.pi
    sum_part = sum_coefficient * pi_value
    
    # Part 2: The Integral
    # As reasoned in the explanation, the complex integral simplifies. Its value is 1.
    integral_part = 1.0
    
    # Final Result
    final_result = sum_part * integral_part
    
    # Print the equation with each number as requested.
    print(f"({sum_coefficient} * {pi_value}) * {integral_part} = {final_result}")

solve_problem()