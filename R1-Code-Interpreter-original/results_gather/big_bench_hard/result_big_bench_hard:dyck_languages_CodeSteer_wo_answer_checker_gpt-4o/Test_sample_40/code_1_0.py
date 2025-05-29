def complete_brackets(sequence):
    # Stack to keep track of opening brackets
    stack = []
    
    # Dictionary to map closing brackets to opening brackets
    bracket_map = {')': '(', '}': '{', '>': '<'}
    
    # Result to store the closing sequence
    closing_sequence = ""
    
    # Iterate through each character in the sequence
    for char in sequence:
        if char in bracket_map.values():
            # If it's an opening bracket, push onto the stack
            stack.append(char)
        elif char in bracket_map:
            # If it's a closing bracket, pop from the stack
            if stack and stack[-1] == bracket_map[char]:
                stack.pop()
    
    # Close any remaining opening brackets in the stack
    closing_map = {v: k for k, v in bracket_map.items()}
    while stack:
        closing_sequence += closing_map[stack.pop()]
    
    # Combine the input sequence with the closing sequence
    complete_sequence = sequence + closing_sequence
    
    # Output the final answer
    return f"<<<{complete_sequence}>>>"

# Input sequence
input_sequence = "{ < { } > { ("

# Get the complete sequence
result = complete_brackets(input_sequence)

# Print the result
print(result)