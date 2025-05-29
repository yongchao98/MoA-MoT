def find_missing_brackets(sequence):
    # Initialize stacks for different types of brackets
    stack = []
    input_str = sequence
    result = sequence
    
    # Dictionary to map closing brackets to opening brackets
    bracket_pairs = {
        '}': '{',
        ']': '[',
        ')': '(',
        '>': '<'
    }
    
    # Dictionary for opposite pairs
    opposite_pairs = {
        '{': '}',
        '[': ']',
        '(': ')',
        '<': '>'
    }
    
    # First pass: check existing sequence and build stack
    for char in input_str:
        if char in '[{(<':
            stack.append(char)
        elif char in ']})>':
            if not stack:
                continue
            if stack[-1] == bracket_pairs[char]:
                stack.pop()
            else:
                continue
    
    # Add closing brackets for remaining opening brackets
    missing_closing = ''
    while stack:
        missing_closing = opposite_pairs[stack.pop()] + missing_closing
    
    print(f"Original sequence: {sequence}")
    print(f"Completed sequence: {sequence + missing_closing}")

# Test with the given sequence
sequence = "[ < [ ] ( ) ( ( { { } } ) ) < { < > } > [ ] > ] < ( ) > ( ( ( ) ) ) ( < >"
find_missing_brackets(sequence)