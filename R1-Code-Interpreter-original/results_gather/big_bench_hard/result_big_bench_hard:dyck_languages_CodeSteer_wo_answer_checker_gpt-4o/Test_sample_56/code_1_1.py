def complete_sequence(sequence):
    # Define matching pairs
    matching_brackets = {'<': '>', '[': ']', '(': ')'}
    
    # Initialize a stack
    stack = []
    
    # Iterate through the sequence
    for char in sequence:
        if char in matching_brackets:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_brackets.values():  # If it's a closing bracket
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()
            else:
                # If there's a mismatch, the sequence is invalid
                return "Invalid sequence"
    
    # If the stack is not empty, add the corresponding closing brackets
    while stack:
        sequence += matching_brackets[stack.pop()]
    
    return f"<<<{sequence}>>>"

# Input sequence
input_sequence = "< < [ ( ) ] >"
# Output the completed sequence
print(complete_sequence(input_sequence))