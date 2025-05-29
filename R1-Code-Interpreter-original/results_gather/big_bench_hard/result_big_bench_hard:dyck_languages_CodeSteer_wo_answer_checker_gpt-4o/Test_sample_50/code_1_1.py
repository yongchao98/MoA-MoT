def complete_brackets(sequence):
    # Define matching pairs for brackets
    matching_bracket = {
        '(': ')',
        '<': '>',
        '[': ']',
        '{': '}'
    }
    
    # Initialize a stack to keep track of opening brackets
    stack = []
    
    # Result list to build the completed sequence
    result = []
    
    # Iterate through each character in the sequence
    for char in sequence:
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
            result.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack:
                stack.pop()  # Pop the matching opening bracket
            result.append(char)
    
    # Close any remaining opening brackets in the stack
    while stack:
        open_bracket = stack.pop()
        result.append(matching_bracket[open_bracket])
    
    # Join the result list into a string
    completed_sequence = ' '.join(result)
    
    # Output the final answer
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "( < > ( [ ( )"

# Complete the sequence
complete_brackets(input_sequence)