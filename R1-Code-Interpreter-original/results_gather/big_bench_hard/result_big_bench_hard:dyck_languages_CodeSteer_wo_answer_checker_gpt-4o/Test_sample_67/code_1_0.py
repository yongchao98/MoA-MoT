def complete_sequence(input_sequence):
    # Define matching pairs of brackets
    matching_brackets = {
        '[': ']',
        '<': '>',
        '(': ')',
        '{': '}'
    }
    
    # Stack to keep track of open brackets
    stack = []
    
    # Process the input sequence
    for char in input_sequence:
        if char in matching_brackets:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_brackets.values():  # If it's a closing bracket
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by closing remaining open brackets
    completed_sequence = input_sequence + ''.join(matching_brackets[char] for char in reversed(stack))
    
    # Output the completed sequence
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "[ [ < [ ( [ < [ ] > ] ) ] > ] ( { }"
complete_sequence(input_sequence)