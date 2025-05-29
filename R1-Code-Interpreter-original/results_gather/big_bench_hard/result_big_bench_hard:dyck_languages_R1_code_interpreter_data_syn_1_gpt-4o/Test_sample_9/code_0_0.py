def complete_sequence(sequence):
    # Define matching pairs
    matching_pairs = {
        ']': '[',
        ')': '(',
        '}': '{',
        '>': '<'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Result to store the completed sequence
    result = []
    
    # Process each character in the sequence
    for char in sequence:
        if char in '[({<':
            # If it's an opening bracket, push to stack
            stack.append(char)
            result.append(char)
        elif char in '])}>':
            # If it's a closing bracket, check if it matches the top of the stack
            if stack and stack[-1] == matching_pairs[char]:
                stack.pop()
                result.append(char)
            else:
                # If it doesn't match, the sequence is invalid
                return "Invalid sequence"
        else:
            # Add other characters to the result
            result.append(char)
    
    # Add the necessary closing brackets for any remaining opening brackets in the stack
    while stack:
        open_bracket = stack.pop()
        for close_bracket, open_bracket_match in matching_pairs.items():
            if open_bracket == open_bracket_match:
                result.append(close_bracket)
                break
    
    return ''.join(result)

# Input sequence
input_sequence = "[ ] ( [ [ { < { { ( < > [ ] ) } } < > > } ] ] { }"

# Complete the sequence
completed_sequence = complete_sequence(input_sequence)

# Print the completed sequence
print(completed_sequence)