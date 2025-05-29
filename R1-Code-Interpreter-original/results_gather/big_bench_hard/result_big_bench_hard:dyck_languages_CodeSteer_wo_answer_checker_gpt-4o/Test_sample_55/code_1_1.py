def complete_sequence(sequence):
    # Count the number of opening and closing symbols
    open_paren = sequence.count('(')
    close_paren = sequence.count(')')
    open_angle = sequence.count('<')
    close_angle = sequence.count('>')

    # Calculate the number of closing symbols needed
    needed_close_paren = open_paren - close_paren
    needed_close_angle = open_angle - close_angle

    # Append the necessary closing symbols
    completed_sequence = sequence + ')' * needed_close_paren + '>' * needed_close_angle

    # Output the completed sequence in the specified format
    print(completed_sequence)

# Input sequence
input_sequence = "( < < < >"
complete_sequence(input_sequence)