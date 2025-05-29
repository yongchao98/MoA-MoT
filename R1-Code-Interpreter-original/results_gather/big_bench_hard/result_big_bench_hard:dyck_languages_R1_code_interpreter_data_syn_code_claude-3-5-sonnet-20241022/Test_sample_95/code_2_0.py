def analyze_and_complete_sequence(s):
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    stack = []
    sequence = s.replace(" ", "")
    result = ""
    
    # First pass: validate and build the sequence
    for char in sequence:
        result += char
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                print(f"Invalid sequence - extra closing bracket: {char}")
                return
            if pairs[stack[-1]] == char:
                stack.pop()
            else:
                print(f"Invalid sequence - mismatched brackets")
                return
    
    # Add remaining closing brackets in reverse order
    while stack:
        result += pairs[stack.pop()]
    
    # Validate final sequence
    validation_stack = []
    for char in result:
        if char in '([{<':
            validation_stack.append(char)
        elif char in ')]}>':
            if not validation_stack or pairs[validation_stack.pop()] != char:
                print("Validation failed")
                return
    
    if not validation_stack:  # If stack is empty, sequence is valid
        print(result)

# Test the input
sequence = "( ( ( [ < [ < [ ( [ { < ( < ( ( ) ) { } [ ] > ) > } ] ) { < > } ] [ < [ { ( < < { { [ < { [ [ < > [ ] ] ] } > ] { { } } } } > > ) } ] > { < > { } } ] > ( [ ] ) ] > ] ) ( )"
analyze_and_complete_sequence(sequence)