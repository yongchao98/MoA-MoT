def complete_sequence(sequence):
    # Define matching pairs
    matching_bracket = {'[': ']', '<': '>', '{': '}', '(': ')'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through the sequence
    for char in sequence:
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
            else:
                # This should not happen as per the problem statement
                return "Invalid sequence"
    
    # Complete the sequence by adding the necessary closing brackets
    completed_sequence = sequence + ''.join(matching_bracket[char] for char in reversed(stack))
    
    # Output the completed sequence
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "[ [ < < { [ ] ( [ [ ( [ [ { [ ] } ] ] ) ( ) ] ] ) ( ) [ ] < [ ( ( ( ) ) ) ] < < < > [ ( ) ] > [ ] > > [ < > ] } > >"

# Call the function with the input sequence
complete_sequence(input_sequence)