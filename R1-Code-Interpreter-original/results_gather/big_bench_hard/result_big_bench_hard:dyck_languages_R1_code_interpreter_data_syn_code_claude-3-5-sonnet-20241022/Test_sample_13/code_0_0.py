def find_missing_closing(sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Mapping of closing to opening brackets
    pairs = {')': '(', '}': '{', ']': '[', '>': '<'}
    # Mapping of opening to closing brackets
    opening_to_closing = {'(': ')', '{': '}', '[': ']', '<': '>'}
    
    # Process the input sequence
    for char in sequence:
        if char in '({[<':
            stack.append(char)
        elif char in ')}]>':
            if not stack or stack[-1] != pairs[char]:
                return "Invalid sequence"
            stack.pop()
    
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing = opening_to_closing[stack.pop()] + missing
    
    print(f"Original sequence: {sequence}")
    print(f"Missing closing brackets: {missing}")
    return missing

# Test the sequence
sequence = "< ( ( [ < > { [ { ( ) } ] < { < { } > [ ( < > ) ] } > } [ < > ] ] ) { { ( { ( ( [ ( [ ] ) ] < { } > ) ) { { ( [ [ ] ] ) } [ ( ) ] { { [ ] } } } } ) } ( { } ) }"
result = find_missing_closing(sequence)