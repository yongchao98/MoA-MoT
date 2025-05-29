def complete_sequence(sequence):
    # Define matching pairs
    matching_brackets = {'(': ')', '{': '}', '<': '>'}
    # Stack to keep track of opening brackets
    stack = []
    # Result list to build the completed sequence
    result = []

    for char in sequence:
        if char in matching_brackets:
            # If it's an opening bracket, push to stack
            stack.append(char)
            result.append(char)
        elif char in matching_brackets.values():
            # If it's a closing bracket, check the stack
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()
                result.append(char)
            else:
                # If there's a mismatch or no opening bracket, ignore the closing bracket
                continue
        else:
            # Add any non-bracket characters directly to the result
            result.append(char)

    # Add the necessary closing brackets for any remaining opening brackets
    while stack:
        result.append(matching_brackets[stack.pop()])

    return ''.join(result)

input_sequence = "( { } ) < { < { } >"
completed_sequence = complete_sequence(input_sequence)
print(f"<<<{completed_sequence}>>>")