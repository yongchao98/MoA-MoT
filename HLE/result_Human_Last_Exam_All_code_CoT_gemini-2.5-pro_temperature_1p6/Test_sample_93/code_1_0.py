import math

def solve_sum_of_squares():
    """
    This function calculates the number of non-negative integer solutions to the Diophantine equation:
    x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024
    """
    
    # The target sum of the equation
    target = 2024
    
    # The number of variables in the sum
    num_vars = 5
    
    # dp[n] will store the number of ways to write n as a sum of k squares.
    # We update this list for k = 1, 2, ..., num_vars.
    dp = [0] * (target + 1)
    
    # Base case: k = 1.
    # The number of solutions to x_1^2 = n is 1 if n is a perfect square, and 0 otherwise.
    limit = int(math.sqrt(target))
    for i in range(limit + 1):
        dp[i * i] = 1
        
    # Iterate for k from 2 to num_vars, building up the solution
    for k in range(2, num_vars + 1):
        # Create a new list to store the results for the sum with k variables
        new_dp = [0] * (target + 1)
        for n in range(target + 1):
            # Calculate dp[k][n] using the recurrence relation:
            # C(k, n) = sum_{i=0 to sqrt(n)} C(k-1, n - i*i)
            limit_i = int(math.sqrt(n))
            for i in range(limit_i + 1):
                # Summing up the solutions for k-1 variables for the remainder n - i^2
                new_dp[n] += dp[n - i * i]
        # Update dp to the newly calculated values for the next iteration
        dp = new_dp
        
    # The final result is the number of solutions for 5 variables to sum to 2024
    result = dp[target]
    
    # Print the result with the numbers from the final equation, as requested.
    print(f"The number of non-negative integer solutions for x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024 is: {result}")

if __name__ == '__main__':
    solve_sum_of_squares()