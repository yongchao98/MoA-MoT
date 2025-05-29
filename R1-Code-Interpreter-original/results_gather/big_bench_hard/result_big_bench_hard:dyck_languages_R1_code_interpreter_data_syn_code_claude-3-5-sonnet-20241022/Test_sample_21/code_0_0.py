def find_missing_brackets(s):
    stack = []
    brackets = {
        '{': '}',
        '(': ')',
        '[': ']'
    }
    
    # Process existing sequence
    for char in s:
        if char in '{([':
            stack.append(char)
        elif char in '})]':
            if not stack:
                return "Invalid sequence"
            if char != brackets[stack.pop()]:
                return "Invalid sequence"
    
    # Find missing closing brackets
    missing = ''
    while stack:
        missing = brackets[stack.pop()] + missing
    
    # Print original sequence + missing brackets
    print(f"Original: {s}")
    print(f"Missing brackets: {missing}")
    print(f"Complete sequence: {s + missing}")

# Test the sequence
find_missing_brackets("{ ( [ ] ) ")