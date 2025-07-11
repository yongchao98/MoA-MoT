def solve_sequence():
    """
    Calculates the next number in the sequence based on the hypothesis
    that it is the difference between the last two numbers.
    """
    sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]

    # The last two numbers in the sequence
    last_number = sequence[-1]
    second_to_last_number = sequence[-2]

    # The proposed equation is to find the difference
    next_number = last_number - second_to_last_number

    # Print out each number in the final equation, as requested
    print(f"{last_number} - {second_to_last_number} = {next_number}")

solve_sequence()