def find_missing_brackets(sequence):
    # Define matching brackets
    brackets = {
        '<': '>',
        '{': '}',
        '[': ']',
        '(': ')'
    }
    
    # Remove spaces and split into characters
    chars = [c for c in sequence if c.strip()]
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in chars:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if stack and brackets[stack[-1]] == char:
                stack.pop()
                
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing += brackets[stack.pop()]
    
    # Print original sequence + missing brackets
    print(f"Original sequence: {sequence}")
    print(f"Missing brackets: {missing}")
    print(f"Complete sequence: {sequence + missing}")

# Test with the given sequence
sequence = "< { < [ [ ( { } ) ] ]"
find_missing_brackets(sequence)