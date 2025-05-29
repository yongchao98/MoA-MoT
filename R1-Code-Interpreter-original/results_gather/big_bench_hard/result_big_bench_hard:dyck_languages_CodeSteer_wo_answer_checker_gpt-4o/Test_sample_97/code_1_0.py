def complete_sequence(sequence):
    # Initialize counters for each type of bracket
    open_round = 0
    close_round = 0
    open_angle = 0
    close_angle = 0

    # Count the brackets in the sequence
    for char in sequence:
        if char == '(':
            open_round += 1
        elif char == ')':
            close_round += 1
        elif char == '<':
            open_angle += 1
        elif char == '>':
            close_angle += 1

    # Calculate the number of closing brackets needed
    needed_close_round = open_round - close_round
    needed_close_angle = open_angle - close_angle

    # Construct the completed sequence
    completed_sequence = sequence + ')' * needed_close_round + '>' * needed_close_angle

    # Output the completed sequence
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "( < >"
complete_sequence(input_sequence)