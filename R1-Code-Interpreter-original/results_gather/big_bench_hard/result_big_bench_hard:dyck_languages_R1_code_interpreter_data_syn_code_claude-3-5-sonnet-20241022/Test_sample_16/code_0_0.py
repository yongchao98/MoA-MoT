def find_missing_brackets(sequence):
    # Define matching pairs
    pairs = {
        '[': ']',
        '<': '>',
        '{': '}'
    }
    
    # Convert string to list of characters
    input_seq = list(sequence.replace(" ", ""))
    
    # Find what brackets are already there
    stack = []
    for char in input_seq:
        if char in '[{<':
            stack.append(char)
        elif char in ']}>':
            if not stack:
                return "Invalid sequence"
            if char == pairs[stack[-1]]:
                stack.pop()
    
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing += pairs[stack.pop()]
    
    # Combine original and missing
    result = ''.join(input_seq) + missing
    
    print(f"Original sequence: {sequence}")
    print(f"Missing brackets: {missing}")
    print(f"Complete sequence: {result}")

# Test the sequence
find_missing_brackets("[[<<{}>")