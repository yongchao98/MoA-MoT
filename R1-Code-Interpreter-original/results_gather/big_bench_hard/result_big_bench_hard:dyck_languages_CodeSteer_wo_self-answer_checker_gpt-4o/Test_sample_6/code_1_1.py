def complete_parentheses(sequence):
    # Define matching pairs
    matching_parentheses = {')': '(', '}': '{', ']': '['}
    opening_parentheses = set(matching_parentheses.values())
    
    # Initialize a stack
    stack = []
    
    # Iterate through the sequence
    for char in sequence:
        if char in opening_parentheses:
            # Push opening parentheses onto the stack
            stack.append(char)
        elif char in matching_parentheses:
            # Check for matching opening parenthesis
            if stack and stack[-1] == matching_parentheses[char]:
                stack.pop()
            else:
                # Invalid sequence
                return "Invalid sequence"
    
    # Construct the completed sequence
    completed_sequence = sequence
    while stack:
        # Pop from stack and find the corresponding closing parenthesis
        opening = stack.pop()
        for closing, opening_match in matching_parentheses.items():
            if opening_match == opening:
                completed_sequence += closing
                break
    
    return f"<<<{completed_sequence}>>>"

# Input sequence
input_sequence = "( { { } }"

# Output the completed sequence
print(complete_parentheses(input_sequence))