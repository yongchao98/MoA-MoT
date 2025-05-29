def complete_sequence(input_sequence):
    # Define matching pairs
    matching_brackets = {'<': '>', '[': ']', '{': '}', '(': ')'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through the input sequence
    for char in input_sequence:
        if char in matching_brackets:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_brackets.values():  # If it's a closing bracket
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()  # Pop the matching opening bracket
    
    # Complete the sequence by adding the necessary closing brackets
    completed_sequence = input_sequence
    while stack:
        completed_sequence += matching_brackets[stack.pop()]
    
    # Output the completed sequence
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "< [ { { < ( ) > { < { } > ( < ( ) > { < [ ( { { ( < [ ] > ) } } { ( ( [ [ { } [ ] ] ] ) ) } ) ] > } ) } }"
complete_sequence(input_sequence)