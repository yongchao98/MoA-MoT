def find_missing_brackets(sequence):
    # Define matching brackets
    brackets = {
        '{': '}',
        '[': ']',
        '(': ')',
        '<': '>'
    }
    
    # Stack to keep track of open brackets
    stack = []
    
    # Process the input sequence
    for char in sequence:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if stack and brackets[stack[-1]] == char:
                stack.pop()
                
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing += brackets[stack.pop()]
        
    # Print original sequence and the missing brackets
    print(f"Original sequence: {sequence}")
    print(f"Missing brackets: {missing}")
    print(f"Complete sequence: {sequence + missing}")

# Test with the given sequence
sequence = "{ { [ { < { } > } ( ) ]"
find_missing_brackets(sequence)