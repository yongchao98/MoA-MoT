def complete_brackets(sequence):
    # Define matching pairs
    pairs = {
        '<': '>', 
        '{': '}',
        '[': ']',
        '(': ')'
    }
    
    # Split the sequence into characters
    chars = [c for c in sequence.strip() if c != ' ']
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in chars:
        if char in pairs.keys():  # Opening bracket
            stack.append(char)
        elif char in pairs.values():  # Closing bracket
            if not stack:
                return "Invalid sequence"
            if char != pairs[stack.pop()]:
                return "Invalid sequence"
    
    # Generate needed closing brackets
    completion = ''
    while stack:
        completion += pairs[stack.pop()]
    
    # Print original + completion
    print(f"Original: {''.join(chars)}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {''.join(chars) + completion}")

# Test the input
sequence = "< { } [ [ ( ( < { < { [ ] } > } > ) )"
complete_brackets(sequence)