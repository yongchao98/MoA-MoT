def complete_brackets(sequence):
    # Define matching brackets
    brackets = {
        '{': '}',
        '(': ')',
        '<': '>',
        '[': ']'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in sequence:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if stack and brackets[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate completion
    completion = ''
    while stack:
        completion += brackets[stack.pop()]
    
    # Print original + completion
    print(f"Original: {sequence}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {sequence + completion}")

# Test the sequence
sequence = "{ ( < [ < > ]"
complete_brackets(sequence)