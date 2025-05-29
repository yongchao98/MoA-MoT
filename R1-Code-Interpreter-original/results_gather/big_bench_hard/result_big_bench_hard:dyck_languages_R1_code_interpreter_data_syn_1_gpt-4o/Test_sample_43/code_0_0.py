def complete_bracket_sequence(sequence):
    open_count = sequence.count('{')
    close_count = sequence.count('}')
    # Calculate the number of closing brackets needed
    needed_closing = open_count - close_count
    # Complete the sequence
    completed_sequence = sequence + '}' * needed_closing
    return completed_sequence

# Given input
input_sequence = "{ { { } }"
# Get the completed sequence
completed_sequence = complete_bracket_sequence(input_sequence)
print(completed_sequence)