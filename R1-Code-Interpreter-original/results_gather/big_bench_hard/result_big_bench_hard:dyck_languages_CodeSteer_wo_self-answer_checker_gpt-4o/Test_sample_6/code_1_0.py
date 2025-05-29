def complete_sequence(sequence):
    # Count the number of each type of parenthesis
    open_paren = sequence.count('(')
    close_paren = sequence.count(')')
    open_brace = sequence.count('{')
    close_brace = sequence.count('}')

    # Calculate the needed closing parentheses
    needed_close_paren = open_paren - close_paren
    needed_close_brace = open_brace - close_brace

    # Complete the sequence
    completed_sequence = sequence + ')' * needed_close_paren + '}' * needed_close_brace

    return completed_sequence

# Input sequence
input_sequence = "( { { } }"

# Complete the sequence
completed_sequence = complete_sequence(input_sequence)

# Output the final answer
print(f"<<<{completed_sequence}>>>")