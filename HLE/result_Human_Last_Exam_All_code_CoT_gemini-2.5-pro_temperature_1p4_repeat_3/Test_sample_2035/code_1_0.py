def solve_problem():
    """
    Calculates the final value based on the derived formula.
    """
    n = 4048
    
    # The term to calculate is sum(alpha_i^2) for i from 1 to n.
    # Given alpha_i = sqrt(i), alpha_i^2 = i.
    # So we need to calculate the sum of the first n integers.
    # The formula for the sum of the first n integers is n * (n + 1) / 2.
    sum_alpha_sq = n * (n + 1) // 2
    
    # The second term in the equation.
    constant_term = 10**15
    
    # The final result is the sum of these two numbers.
    final_result = sum_alpha_sq + constant_term
    
    # Print the numbers in the final equation as requested.
    print(f"The equation is: {sum_alpha_sq} + {constant_term}")
    print(f"The final result is: {final_result}")

solve_problem()