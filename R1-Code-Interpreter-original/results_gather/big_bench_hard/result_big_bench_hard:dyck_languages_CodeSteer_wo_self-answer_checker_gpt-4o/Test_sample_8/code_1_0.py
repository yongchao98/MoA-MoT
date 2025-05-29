def complete_sequence(sequence):
    stack = []
    completed_sequence = []

    # Mapping of opening to closing brackets
    matching_bracket = {'[': ']', '<': '>'}

    for char in sequence:
        completed_sequence.append(char)
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
            else:
                # Unbalanced sequence, but we assume input is correct up to this point
                pass

    # Add necessary closing brackets for any unmatched opening brackets
    while stack:
        completed_sequence.append(matching_bracket[stack.pop()])

    return ''.join(completed_sequence)

# Input sequence
input_sequence = "[ < > ] [ [ < > "

# Complete the sequence
completed_sequence = complete_sequence(input_sequence)

# Output the final answer in the specified format
print(f"<<<{completed_sequence}>>>")