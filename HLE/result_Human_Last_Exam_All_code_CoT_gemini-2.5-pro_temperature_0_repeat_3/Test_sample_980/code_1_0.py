def solve_sequence():
    """
    This function calculates the next number in the sequence based on the
    hypothesis that it is the sum of the last two numbers.
    """
    # The given sequence
    sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]

    # The last two numbers in the sequence
    second_to_last_number = sequence[-2]
    last_number = sequence[-1]

    # The next number is the sum of the last two
    next_number = second_to_last_number + last_number

    # Print the equation as requested
    print(f"The final equation is: {second_to_last_number} + {last_number} = {next_number}")
    print(f"The next number in the sequence is: {next_number}")

solve_sequence()