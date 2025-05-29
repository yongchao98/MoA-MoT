def complete_sequence(sequence):
    # Define matching pairs
    matching_brackets = {'<': '>', '[': ']', '{': '}', '(': ')'}
    
    # Initialize a stack to keep track of opening brackets
    stack = []
    
    # Iterate through each character in the sequence
    for char in sequence:
        if char in matching_brackets:
            # If it's an opening bracket, push it onto the stack
            stack.append(char)
        elif char in matching_brackets.values():
            # If it's a closing bracket, check if it matches the top of the stack
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()
            else:
                # If it doesn't match, the sequence is already invalid
                return "Invalid sequence"
    
    # Construct the completed sequence
    completed_sequence = sequence
    while stack:
        # Append the corresponding closing bracket for each remaining opening bracket
        completed_sequence += matching_brackets[stack.pop()]
    
    return f"<<<{completed_sequence}>>>"

# Input sequence
input_sequence = "< [ ] { < ( ) > } [ ] ( { }"

# Output the completed sequence
print(complete_sequence(input_sequence))