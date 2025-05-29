def complete_sequence(input_sequence):
    # Define matching pairs
    matching_brackets = {'{': '}', '(': ')', '[': ']', '<': '>'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Traverse the input sequence
    for char in input_sequence:
        if char in matching_brackets:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_brackets.values():  # If it's a closing bracket
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Construct the final sequence by adding necessary closing brackets
    final_sequence = input_sequence
    while stack:
        final_sequence += matching_brackets[stack.pop()]
    
    return f"<<<{final_sequence}>>>"

# Input sequence
input_sequence = "{ { ( { [ [ { < [ ] > } ] ] } ) [ ( < { < > } > ) ]"

# Get the completed sequence
completed_sequence = complete_sequence(input_sequence)

# Print the result
print(completed_sequence)