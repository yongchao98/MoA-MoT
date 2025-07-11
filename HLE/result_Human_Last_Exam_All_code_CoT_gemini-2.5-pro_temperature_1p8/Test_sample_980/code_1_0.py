def find_next_in_sequence():
    """
    Calculates the next number in the sequence based on the hypothesis
    that it is the XOR of the last two numbers.
    """
    sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]

    # The last two numbers in the sequence are used for the calculation.
    # The last number is sequence[-1] and the second to last is sequence[-2].
    second_to_last_num = sequence[-2]
    last_num = sequence[-1]

    # The hypothesized rule is: next_number = second_to_last_number XOR last_number
    next_num = second_to_last_num ^ last_num

    print("The proposed final equation is based on the last two numbers of the sequence.")
    print(f"{second_to_last_num} ^ {last_num} = {next_num}")

find_next_in_sequence()