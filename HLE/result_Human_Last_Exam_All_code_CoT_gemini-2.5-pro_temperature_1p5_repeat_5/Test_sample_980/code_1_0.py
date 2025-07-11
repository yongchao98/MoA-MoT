def find_next_number_in_sequence():
    """
    This function analyzes the given sequence to find the next number.
    It looks for an alternating pattern of the form 'x, y, x' at the end of the sequence
    and predicts 'y' as the next number.
    """
    sequence = [111, 142, 111, 41, 67, 67, 67, 93, 111, 111, 62, 62, 111, 111, 36, 36, 49, 155, 49, 62, 49, 49, 62, 62, 10, 36, 36, 36, 124, 124, 124, 36, 124]

    # Check if the sequence has at least 3 elements
    if len(sequence) >= 3:
        # Get the last three numbers
        last_three = sequence[-3:]
        x = last_three[0]
        y = last_three[1]
        
        # Check if the pattern is x, y, x
        if last_three[2] == x and x != y:
            next_number = y
            print("The pattern at the end of the sequence is an alternation.")
            # The prompt asks to "output each number in the final equation"
            # We interpret this as showing the logic of the pattern.
            print(f"The detected pattern is: {x}, {y}, {x}, ...")
            print(f"Based on this alternating pattern, the next number should be {next_number}.")
        else:
            print("Could not determine a clear alternating pattern at the end.")
    else:
        print("Sequence is too short to find a pattern.")

find_next_number_in_sequence()
