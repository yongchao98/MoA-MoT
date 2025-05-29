def complete_sequence(sequence):
    # Define matching pairs
    matching_brackets = {'(': ')', '{': '}', '[': ']', '<': '>'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through the sequence
    for char in sequence:
        if char in matching_brackets:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_brackets.values():  # If it's a closing bracket
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Complete the sequence by adding the necessary closing brackets
    completed_sequence = sequence
    while stack:
        completed_sequence += matching_brackets[stack.pop()]
    
    return f"<<<{completed_sequence}>>>"

# Input sequence
input_sequence = "( < { < { } > { { ( [ ( { } ) ( ( [ ( [ ( { [ { } ] } ) ] ) ] ) ) [ < [ [ [ [ [ < > ] [ { [ [ { ( ( < [ ] > ) [ ( [ ] ) ] < { [ ] } > ) } ] ] } ] { < ( < > ) > } ] ] ] ] > ] ] ) } } } >"

# Output the completed sequence
print(complete_sequence(input_sequence))