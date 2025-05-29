def find_missing_brackets(sequence):
    # Remove spaces
    sequence = sequence.replace(" ", "")
    
    # Define bracket pairs
    brackets = {
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Stack to track opening brackets
    stack = []
    
    # Process existing sequence
    for char in sequence:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if not stack:
                return "Invalid sequence"
            if char != brackets[stack.pop()]:
                return "Invalid sequence"
    
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing += brackets[stack.pop()]
    
    # Print original sequence with completion
    print(f"Original: {sequence}")
    print(f"Completion: {missing}")
    print(f"Full sequence: {sequence + missing}")

# Test the input
find_missing_brackets("[ < [ ] > ] [ { [ ] }")