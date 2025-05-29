def complete_sequence(s):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to match opening and closing brackets
    pairs = {
        '<': '>', 
        '[': ']', 
        '{': '}', 
        '(': ')'
    }
    
    # Process the given sequence
    for char in s:
        if char in '<[{(':
            stack.append(char)
        elif char in '>]}):':
            if not stack:
                return "Invalid sequence"
            if char != pairs[stack.pop()]:
                return "Invalid sequence"
    
    # Generate completion by closing all remaining brackets
    completion = ''
    while stack:
        completion = pairs[stack.pop()] + completion
    
    # Print original sequence and completion
    print(f"Original sequence: {s}")
    print(f"Needed closing brackets: {completion}")
    print(f"Complete sequence: {s + completion}")

# Test the sequence
sequence = "< [ ] { < ( ) > } [ ] ( { }"
complete_sequence(sequence)