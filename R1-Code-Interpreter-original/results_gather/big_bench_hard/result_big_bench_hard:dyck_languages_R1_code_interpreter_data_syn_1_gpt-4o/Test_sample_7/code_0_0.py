def complete_sequence(input_sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to match opening and closing brackets
    matching_bracket = {'(': ')', '[': ']', '{': '}', '<': '>'}
    
    # Iterate over the input sequence
    for char in input_sequence:
        if char in matching_bracket:
            # If it's an opening bracket, push to stack
            stack.append(char)
        elif char in matching_bracket.values():
            # If it's a closing bracket, pop from stack
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
    
    # Generate the closing sequence
    closing_sequence = ''.join(matching_bracket[char] for char in reversed(stack))
    return closing_sequence

# Input sequence
input_sequence = "< ( ( ( [ { } ] )"
# Remove spaces for processing
input_sequence = input_sequence.replace(" ", "")
# Get the closing sequence
closing_sequence = complete_sequence(input_sequence)
print(closing_sequence)