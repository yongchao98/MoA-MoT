def find_missing_brackets(sequence):
    stack = []
    pairs = {'{': '}', '(': ')'}
    
    # Process existing sequence
    for char in sequence:
        if char in '{(':  # Opening brackets
            stack.append(char)
        elif char in '})':  # Closing brackets
            if stack and pairs[stack[-1]] == char:
                stack.pop()
    
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing += pairs[stack.pop()]
    
    # Print original sequence + missing brackets
    print(f"Original sequence: {sequence}")
    print(f"Missing brackets: {missing}")
    print(f"Complete sequence: {sequence + missing}")

# Test with the given sequence
sequence = "{ { } ( ( ) )"
sequence = sequence.replace(" ", "")  # Remove spaces
find_missing_brackets(sequence)