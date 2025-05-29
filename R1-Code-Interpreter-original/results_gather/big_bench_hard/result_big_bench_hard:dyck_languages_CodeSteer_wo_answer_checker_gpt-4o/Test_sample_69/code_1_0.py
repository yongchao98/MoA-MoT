def complete_sequence(input_sequence):
    # Define matching pairs
    matching_bracket = {'{': '}', '<': '>', '[': ']'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Traverse the input sequence
    for char in input_sequence:
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()  # Pop the matching opening bracket
            else:
                return "Invalid sequence"  # Mismatched or unbalanced
    
    # If stack is not empty, we need to add closing brackets
    completion = ''.join(matching_bracket[char] for char in reversed(stack))
    
    # Construct the final answer
    final_answer = input_sequence + completion
    return f"<<<{final_answer}>>>"

# Input sequence
input_sequence = "{ ( { { < { < > } > } < { } > < < < [ < [ ] > ] > > >"

# Get the completed sequence
result = complete_sequence(input_sequence)

# Print the result
print(result)