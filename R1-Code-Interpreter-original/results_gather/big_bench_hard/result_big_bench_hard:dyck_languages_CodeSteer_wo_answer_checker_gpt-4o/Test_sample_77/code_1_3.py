def complete_sequence(sequence):
    stack = []
    completed_sequence = []

    # Mapping of opening to closing brackets
    matching_bracket = {'<': '>', '[': ']', '(': ')'}

    for char in sequence:
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
            completed_sequence.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
            completed_sequence.append(char)
        else:
            completed_sequence.append(char)
        
        # Debug: Print the current state of the stack
        print(f"Current character: {char}, Stack: {stack}")

    # Add the necessary closing brackets for any remaining open brackets
    while stack:
        completed_sequence.append(matching_bracket[stack.pop()])

    return ''.join(completed_sequence)

sequence = "< < [ ] > ( < > ) <"
completed_sequence = complete_sequence(sequence)
print(f"<<<{completed_sequence}>>>")