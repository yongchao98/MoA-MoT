def complete_sequence(input_sequence):
    # Define matching pairs for brackets
    matching_bracket = {'[': ']', '<': '>', '{': '}'}
    
    # Initialize a stack to keep track of opening brackets
    stack = []
    
    # Process the input sequence
    for char in input_sequence:
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()  # Properly matched, pop from stack
    
    # Complete the sequence by adding the necessary closing brackets
    completed_sequence = input_sequence + ''.join(matching_bracket[char] for char in reversed(stack))
    
    return completed_sequence

# Initial sequence
input_sequence = "[ [ < < { } >"

# Get the completed sequence
completed_sequence = complete_sequence(input_sequence)

# Output the final answer
print(f"<<<{completed_sequence}>>>")