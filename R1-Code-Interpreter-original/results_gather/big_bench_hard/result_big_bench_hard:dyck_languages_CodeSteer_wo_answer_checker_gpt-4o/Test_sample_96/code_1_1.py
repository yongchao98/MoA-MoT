def complete_sequence(input_sequence):
    # Define matching pairs of brackets
    matching_brackets = {'<': '>', '{': '}'}
    
    # Initialize a stack to track open brackets
    stack = []
    
    # Process each character in the input sequence
    for char in input_sequence:
        if char in matching_brackets:  # If it's an open bracket
            stack.append(char)
        elif char in matching_brackets.values():  # If it's a close bracket
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()  # Properly closed, pop the stack
            else:
                # Unbalanced, add the corresponding open bracket to the stack
                for open_bracket, close_bracket in matching_brackets.items():
                    if close_bracket == char:
                        stack.append(open_bracket)
                        break
    
    # Construct the completed sequence
    completed_sequence = input_sequence
    while stack:
        completed_sequence += matching_brackets[stack.pop()]
    
    # Output the final answer in the specified format
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "< < { } >"

# Call the function to complete the sequence
complete_sequence(input_sequence)