def solve_sequence():
    """
    This function solves the number sequence puzzle by identifying a two-level pattern.
    1. It establishes a primary relationship: x_n = 3 * x_{n-1} - S_n.
    2. It finds the sequence of subtractions, S_n = [10, 18, 34, ...].
    3. It identifies that the differences in S_n are a geometric progression (8, 16, ...).
    4. It predicts the next subtraction term and calculates the final term in the original sequence.
    """
    # The given sequence
    sequence = [2, 11, 23, 51, 119]

    # The last known term
    last_term = sequence[-1]

    # Calculate the sequence of subtracted numbers (S_n)
    # S_3 = 3 * 11 - 23 = 10
    # S_4 = 3 * 23 - 51 = 18
    # S_5 = 3 * 51 - 119 = 34
    subtractions = [3 * sequence[i-1] - sequence[i] for i in range(2, len(sequence))]

    # The last calculated subtraction
    last_subtraction = subtractions[-1]

    # Calculate the differences in the subtraction sequence
    # diff_1 = 18 - 10 = 8
    # diff_2 = 34 - 18 = 16
    differences = [subtractions[i] - subtractions[i-1] for i in range(1, len(subtractions))]
    
    # The last calculated difference
    last_difference = differences[-1]

    # The pattern in the differences is that they double each time.
    # So the next difference is 16 * 2 = 32.
    next_difference = last_difference * 2

    # The next subtraction number is the last one plus the next difference.
    # 34 + 32 = 66
    next_subtraction = last_subtraction + next_difference

    # The next term in the main sequence is 3 * the last term - the next subtraction number.
    # 3 * 119 - 66
    next_term = 3 * last_term - next_subtraction

    # Print the final equation
    print(f"The pattern is x_n = 3 * x_(n-1) - S_n, where the differences in S_n double.")
    print(f"The sequence of subtractions is {subtractions}, and the next subtraction is {next_subtraction}.")
    print(f"The final calculation is:")
    print(f"3 * {last_term} - {next_subtraction} = {next_term}")

solve_sequence()
<<<A>>>