def complete_sequence(s):
    # Define matching brackets
    brackets = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Process existing sequence
    stack = []
    input_sequence = s.replace(" ", "")
    
    # First validate existing sequence
    for char in input_sequence:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if not stack:  # Stack empty but found closing bracket
                return "Invalid sequence"
            if char != brackets[stack[-1]]:  # Mismatched brackets
                return "Invalid sequence"
            stack.pop()
    
    # Now complete the sequence
    completion = ""
    while stack:
        completion = brackets[stack.pop()] + completion
    
    print(f"Original: {input_sequence}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {input_sequence + completion}")

# Test with the given sequence
test_sequence = "[ < [ { { < > [ ( ) ] ( ( ) ) } } ] > ] { { { ( < > )"
complete_sequence(test_sequence)