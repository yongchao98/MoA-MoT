def find_missing_brackets(sequence):
    # Define matching pairs
    pairs = {
        '[': ']',
        '{': '}',
        '<': '>',
        '(': ')'
    }
    
    # Convert sequence to list for easier manipulation
    sequence = list(sequence.replace(" ", ""))
    
    # Stack to keep track of opening brackets
    stack = []
    
    # First pass: check what brackets are still open
    for char in sequence:
        if char in pairs.keys():  # Opening bracket
            stack.append(char)
        elif char in pairs.values():  # Closing bracket
            if stack and pairs[stack[-1]] == char:
                stack.pop()
    
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing = pairs[stack.pop()] + missing
    
    # Combine original sequence with missing brackets
    result = ''.join(sequence) + missing
    print(f"Complete sequence: {result}")
    print(f"Missing brackets: {missing}")

# Test the sequence
find_missing_brackets("[ < { [ ] } >")