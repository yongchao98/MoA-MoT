def complete_sequence(s):
    # Stack to keep track of opening brackets
    stack = []
    # Original sequence
    sequence = s
    
    # Mapping of opening to closing brackets
    brackets = {
        '<': '>',
        '(': ')',
        '{': '}',
        '[': ']'
    }
    
    # Process existing sequence
    for char in sequence:
        if char in '<({[':
            stack.append(char)
        elif char in '>)}]':
            if not stack:
                return "Invalid sequence"
            if char != brackets[stack.pop()]:
                return "Invalid sequence"
    
    # Complete the sequence by adding closing brackets
    completion = ''
    while stack:
        completion += brackets[stack.pop()]
    
    # Print original + completion
    print(f"Original: {sequence}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {sequence + completion}")

# Test with the given sequence
complete_sequence("< ( ( { < > } )")