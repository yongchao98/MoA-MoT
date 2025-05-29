def complete_bracket_sequence(sequence):
    # Define matching pairs of brackets
    matching_bracket = {'{': '}', '(': ')', '<': '>', '[': ']'}
    
    # Initialize a stack to keep track of opening brackets
    stack = []
    
    # Process each character in the input sequence
    for char in sequence:
        if char in matching_bracket:
            # If it's an opening bracket, push it onto the stack
            stack.append(char)
        elif char in matching_bracket.values():
            # If it's a closing bracket, pop from the stack
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by closing remaining open brackets
    completed_sequence = sequence
    while stack:
        completed_sequence += matching_bracket[stack.pop()]
    
    return completed_sequence

# Input sequence
input_sequence = "{ ( < [ < > ]"

# Get the completed sequence
completed_sequence = complete_bracket_sequence(input_sequence)

# Output the final answer
print(f"<<<{completed_sequence}>>>")