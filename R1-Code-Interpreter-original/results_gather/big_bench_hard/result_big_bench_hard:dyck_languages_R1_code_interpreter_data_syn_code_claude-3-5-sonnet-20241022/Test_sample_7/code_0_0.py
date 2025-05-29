def find_missing_brackets(sequence):
    # Define matching pairs
    pairs = {
        '<': '>',
        '(': ')',
        '[': ']',
        '{': '}'
    }
    
    # Remove spaces and convert to list
    sequence = list(sequence.replace(" ", ""))
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in sequence:
        if char in pairs:  # Opening bracket
            stack.append(char)
        elif char in pairs.values():  # Closing bracket
            if not stack:
                return "Invalid sequence"
            if char != pairs[stack[-1]]:
                return "Invalid sequence"
            stack.pop()
    
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing += pairs[stack.pop()]
    
    # Combine original sequence with missing brackets
    result = ''.join(sequence) + missing
    print(f"Complete sequence: {result}")
    return result

# Test the sequence
sequence = "< ( ( ( [ { } ] )"
find_missing_brackets(sequence)