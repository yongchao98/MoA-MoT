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
    for char in s.replace(" ", ""):
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                return "Invalid sequence"
            if pairs[stack[-1]] != char:
                return "Invalid sequence"
            stack.pop()
    
    # Complete the sequence by adding closing brackets
    completion = ''
    while stack:
        completion += pairs[stack.pop()]
    
    # Print original + completion
    print(s.strip() + completion)

# Test the input
input_seq = "( ( ( [ < [ < [ ( [ { < ( < ( ( ) ) { } [ ] > ) > } ] ) { < > } ] [ < [ { ( < < { { [ < { [ [ < > [ ] ] ] } > ] { { } } } } > > ) } ] > { < > { } } ] > ( [ ] ) ] > ] ) ( )"
complete_sequence(input_seq)