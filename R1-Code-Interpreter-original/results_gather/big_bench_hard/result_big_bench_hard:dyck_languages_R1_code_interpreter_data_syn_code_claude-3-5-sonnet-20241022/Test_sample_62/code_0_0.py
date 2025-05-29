def complete_sequence(s):
    # Define matching pairs
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in s:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>' and stack:
            if pairs[stack[-1]] == char:
                stack.pop()
    
    # Generate completion
    completion = ''
    while stack:
        completion = pairs[stack.pop()] + completion
    
    # Print original sequence and completion
    print(f"Original sequence: {s}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {s + completion}")

# Test with the given input
input_seq = "( ( < [ { { } } ] > { } ) ( < > ) [ ( < > ) ] <"
complete_sequence(input_seq)