def complete_sequence(input_sequence):
    # Define matching pairs
    matching_bracket = {'(': ')', '[': ']', '{': '}', '<': '>'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through the input sequence
    for char in input_sequence:
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
            else:
                # If there's a mismatch, the input is invalid
                return "Invalid input sequence"
    
    # Construct the closing sequence
    closing_sequence = ''.join(matching_bracket[char] for char in reversed(stack))
    
    # Output the final answer
    print(f"<<<{input_sequence}{closing_sequence}>>>")

# Input sequence
input_sequence = "(< < > ( [ < { { } } ( ) > ] ) [ ( ) ] [ ] [ { } ] >"

# Call the function
complete_sequence(input_sequence)