def complete_sequence(sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to match opening and closing brackets
    matching_bracket = {'{': '}', '(': ')', '<': '>', '[': ']'}
    
    # Iterate over the sequence
    for char in sequence:
        if char in matching_bracket:
            # If it's an opening bracket, push to stack
            stack.append(char)
        elif char in matching_bracket.values():
            # If it's a closing bracket, pop from stack
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
    
    # Now, stack contains unclosed opening brackets
    # We need to close them in reverse order
    result = sequence
    while stack:
        result += matching_bracket[stack.pop()]
    
    return result

# Given input sequence
input_sequence = "{ ( < [ < > ]"
# Complete the sequence
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)