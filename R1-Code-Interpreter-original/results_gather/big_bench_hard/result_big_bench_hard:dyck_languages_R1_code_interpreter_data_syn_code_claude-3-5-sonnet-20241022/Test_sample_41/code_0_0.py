def find_completion(sequence):
    # Define matching pairs
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process the input sequence
    for char in sequence:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                return "Invalid sequence"
            if pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate completion
    completion = ''
    while stack:
        completion = pairs[stack.pop()] + completion
    
    # Print both the input and completion
    print(f"Input sequence: {sequence}")
    print(f"Required completion: {completion}")
    return completion

# Test with the given sequence
sequence = "( [ [ { < > { } { } < < < < { } > > ( < ( ( ) ) [ ( [ ] ) ] > ) > > { } } [ < < > > ] < ( ) > ] ]"
find_completion(sequence)