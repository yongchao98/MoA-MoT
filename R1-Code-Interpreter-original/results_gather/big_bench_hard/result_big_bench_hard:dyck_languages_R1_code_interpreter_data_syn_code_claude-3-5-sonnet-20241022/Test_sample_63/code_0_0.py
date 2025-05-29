def find_completion(s):
    # Define matching pairs
    pairs = {
        '[': ']',
        '{': '}',
        '(': ')',
        '<': '>'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in s:
        if char in '[{(<':
            stack.append(char)
        elif char in ']})>':
            if not stack:
                return "Invalid sequence"
            if char == pairs[stack[-1]]:
                stack.pop()
    
    # Generate completion
    completion = ''
    while stack:
        completion = pairs[stack.pop()] + completion
    
    # Print original sequence and needed completion
    print(f"Original sequence: {s}")
    print(f"Needed completion: {completion}")
    print(f"Complete sequence: {s + completion}")

# Test with the given sequence
sequence = "[ < < { } > < > > ( )"
find_completion(sequence)