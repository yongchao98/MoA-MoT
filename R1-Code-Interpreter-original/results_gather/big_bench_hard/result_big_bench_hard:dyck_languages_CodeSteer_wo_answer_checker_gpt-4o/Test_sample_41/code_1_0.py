def complete_sequence(sequence):
    # Define matching pairs
    matching_pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Traverse the sequence
    for char in sequence:
        if char in matching_pairs:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_pairs.values():  # If it's a closing bracket
            if stack and matching_pairs[stack[-1]] == char:
                stack.pop()
            else:
                # Invalid sequence, but we assume input is valid up to this point
                return "Invalid sequence"
    
    # Add the necessary closing brackets
    completed_sequence = sequence
    while stack:
        completed_sequence += matching_pairs[stack.pop()]
    
    return f"<<<{completed_sequence}>>>"

# Input sequence
input_sequence = "( [ [ { < > { } { } < < < < { } > > ( < ( ( ) ) [ ( [ ] ) ] > ) > > { } } [ < < > > ] < ( ) > ] ]"

# Get the completed sequence
result = complete_sequence(input_sequence)

# Print the result
print(result)