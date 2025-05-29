def complete_sequence(input_sequence):
    # Define matching pairs
    matching_bracket = {
        '}': '{',
        '>': '<',
        ']': '[',
        ')': '('
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through the input sequence
    for char in input_sequence:
        if char in '{[<(':
            # Push opening brackets onto the stack
            stack.append(char)
        elif char in '}])>':
            # Check if the stack is empty or the top of the stack doesn't match
            if not stack or stack[-1] != matching_bracket[char]:
                return "Invalid sequence"
            # Pop the matching opening bracket
            stack.pop()
    
    # If the stack is not empty, we need to add closing brackets
    closing_sequence = ''
    while stack:
        opening_bracket = stack.pop()
        for closing, opening in matching_bracket.items():
            if opening == opening_bracket:
                closing_sequence += closing
                break
    
    return closing_sequence

# Input sequence
input_sequence = "{ { { < { } < ( ( < > ) ) > > } } { [ [ { } ] ] [ ] < < { } > [ ( ) { ( { } ) } ] > ( [ < [ ( ( { } ( ) < > < [ < < { { ( [ { } ] ) } } { } { } > > ] > ( ) ) ) < > ] > ] )"

# Get the output
output = complete_sequence(input_sequence)
print(output)