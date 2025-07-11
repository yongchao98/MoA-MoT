def solve_sequence():
    """
    This function solves for the next number in a given integer sequence
    by analyzing the sequence of differences.
    """
    # The given sequence
    sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]

    # Calculate the differences between consecutive numbers
    differences = [sequence[i+1] - sequence[i] for i in range(len(sequence) - 1)]

    # The problem is to predict the next difference.
    # Let's look at the end of the difference sequence:
    # The last few values are [..., 88, 0, 0, -88, 88]
    # We can observe patterns like (x, -x) or (-x, x) pairs, often followed by 0.
    # e.g., (106, -106) is followed by 13. (13, -13) is followed by 0.
    # The sequence also has many zeros.
    # A plausible prediction for the next difference, following the pair (-88, 88), is 0.
    next_difference = 0

    # The last number in the original sequence
    last_number = sequence[-1]

    # Calculate the next number in the sequence
    next_number = last_number + next_difference
    
    # The final equation to find the next number
    print(f"The original sequence is: {sequence}")
    print(f"The sequence of differences is: {differences}")
    print(f"The last number is {last_number}.")
    print(f"The predicted next difference is {next_difference}.")
    print("The equation for the next number is:")
    print(f"{last_number} + {next_difference} = {next_number}")

solve_sequence()