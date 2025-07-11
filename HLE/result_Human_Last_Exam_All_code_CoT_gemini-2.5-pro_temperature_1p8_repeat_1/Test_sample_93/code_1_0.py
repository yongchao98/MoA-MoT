import collections

def count_sum_of_squares_solutions():
    """
    Calculates the number of non-negative integer solutions to the Diophantine equation:
    x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024
    """
    
    target_sum = 2024
    num_variables = 5
    
    # The maximum value for any x_i is floor(sqrt(2024)) = 44.
    max_val = int(target_sum**0.5)
    squares = [i * i for i in range(max_val + 1)]
    
    # dp[k] will store the number of ways to write k as a sum of i squares.
    # We initialize it for i=1. A sum 's' can be formed in 1 way if s is a perfect square.
    dp = collections.defaultdict(int)
    for s in squares:
        if s <= target_sum:
            dp[s] += 1
            
    # Now, we iteratively compute the counts for sums of 2, 3, 4, and 5 squares.
    for i in range(2, num_variables + 1):
        next_dp = collections.defaultdict(int)
        # The current dp has counts for sums of (i-1) squares.
        # We convolve it with the list of single squares to get counts for i squares.
        for current_sum, count in dp.items():
            for s in squares:
                new_sum = current_sum + s
                if new_sum <= target_sum:
                    next_dp[new_sum] += count
        dp = next_dp
        
    result = dp[target_sum]
    
    # The problem asks to "output each number in the final equation".
    # We will print the equation and the final count.
    # The numbers in the equation are: 5 variables, power of 2, and the sum 2024.
    
    variable_terms = " + ".join([f"x_{i}^2" for i in range(1, num_variables + 1)])
    equation = f"{variable_terms} = {target_sum}"
    
    print(f"For the Diophantine equation: {equation}")
    print(f"The number of non-negative integer solutions is: {result}")

count_sum_of_squares_solutions()
<<<100880>>>