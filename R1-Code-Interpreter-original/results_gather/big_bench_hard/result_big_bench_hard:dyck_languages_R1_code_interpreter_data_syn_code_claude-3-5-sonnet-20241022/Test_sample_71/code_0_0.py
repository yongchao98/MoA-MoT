def find_closing_sequence(input_str):
    # Define matching pairs
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process input string
    for char in input_str:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack or pairs[stack.pop()] != char:
                return "Invalid sequence"
    
    # Generate closing sequence
    closing_sequence = ''
    while stack:
        closing_sequence = pairs[stack.pop()] + closing_sequence
    
    # Print original + completion
    print(f"Original: {input_str}")
    print(f"Completion needed: {closing_sequence}")
    print(f"Complete sequence: {input_str + closing_sequence}")
    return closing_sequence

# Test with the given input
input_str = "( ( < { [ [ { [ [ ( [ < < > > ( ( [ ( ) ( { ( < [ ] > ) } ) ] ) ) ] ) ] ( ) ] } { [ [ ] ] } ] ] } > [ { < > } ]"
result = find_closing_sequence(input_str)