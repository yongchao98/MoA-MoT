def complete_sequence(input_sequence):
    # Define matching pairs
    matching_bracket = {')': '(', '}': '{', ']': '[', '>': '<'}
    opening_brackets = set(matching_bracket.values())
    
    # Initialize a stack
    stack = []
    
    # Iterate over each character in the input sequence
    for char in input_sequence:
        if char in opening_brackets:
            # If it's an opening bracket, push it onto the stack
            stack.append(char)
        elif char in matching_bracket:
            # If it's a closing bracket, check if it matches the top of the stack
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                # If it doesn't match, the sequence is invalid
                return "Invalid sequence"
    
    # If the stack is not empty, it means there are unmatched opening brackets
    if stack:
        return "Invalid sequence"
    
    # If we reach here, the sequence is valid
    return "Valid sequence"

# Input sequence
input_sequence = "[ { < > < [ ] > [ < { } > ] [ ( { } ) ] [ ] [ ( ) ] ( ) < { { < { { { { { ( { } ) } ( { < > } ) } } } } > } }"

# Check if the sequence is valid
result = complete_sequence(input_sequence)
print(result)