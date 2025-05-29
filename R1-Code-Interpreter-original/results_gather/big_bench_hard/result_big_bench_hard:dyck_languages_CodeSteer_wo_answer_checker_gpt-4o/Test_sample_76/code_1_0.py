def complete_sequence(sequence):
    # Define matching pairs
    matching_bracket = {'(': ')', '{': '}', '<': '>'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Traverse the input sequence
    for char in sequence:
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by adding the necessary closing brackets
    completed_sequence = sequence
    while stack:
        completed_sequence += matching_bracket[stack.pop()]
    
    # Output the completed sequence
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "( { } < > { < { ( < > ) } >"
complete_sequence(input_sequence)