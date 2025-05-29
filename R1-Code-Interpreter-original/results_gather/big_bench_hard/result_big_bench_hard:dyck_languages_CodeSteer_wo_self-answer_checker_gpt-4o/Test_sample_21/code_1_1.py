def complete_bracket_sequence(sequence):
    # Dictionary to match opening and closing brackets
    matching_bracket = {')': '(', ']': '[', '}': '{'}
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through each character in the sequence
    for char in sequence:
        if char in matching_bracket.values():
            # If it's an opening bracket, push onto the stack
            stack.append(char)
        elif char in matching_bracket:
            # If it's a closing bracket, check for matching opening bracket
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                # If no match, the sequence is invalid
                return "Invalid sequence"
    
    # Add the necessary closing brackets for any remaining opening brackets
    completed_sequence = sequence + ''.join({v: k for k, v in matching_bracket.items()}[bracket] for bracket in reversed(stack))
    
    # Output the final answer in the required format
    return f"<<<{completed_sequence}>>>"

# Given input sequence
input_sequence = "{ ( [ ] )"

# Complete the sequence
print(complete_bracket_sequence(input_sequence))