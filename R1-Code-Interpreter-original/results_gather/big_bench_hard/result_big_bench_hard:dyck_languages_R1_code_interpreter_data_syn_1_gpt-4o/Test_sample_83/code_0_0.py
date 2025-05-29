def complete_sequence(sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to map opening brackets to their corresponding closing brackets
    bracket_map = {'{': '}', '<': '>', '[': ']', '(': ')'}
    
    # Iterate over each character in the sequence
    for char in sequence:
        if char in bracket_map:
            # If it's an opening bracket, push to stack
            stack.append(char)
        elif char in bracket_map.values():
            # If it's a closing bracket, pop from stack
            if stack and bracket_map[stack[-1]] == char:
                stack.pop()
    
    # Now, complete the sequence by adding the necessary closing brackets
    completed_sequence = sequence
    while stack:
        completed_sequence += bracket_map[stack.pop()]
    
    return completed_sequence

# Given input sequence
input_sequence = "{ < [ ]"
# Complete the sequence
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)