def find_missing_brackets(sequence):
    # Define matching pairs
    pairs = {
        '[': ']',
        '{': '}',
        '(': ')',
        '<': '>'
    }
    
    # Convert input to list for processing
    input_seq = [char for char in sequence if char != ' ']
    
    # Stack to track opening brackets
    stack = []
    
    # Process existing sequence
    for char in input_seq:
        if char in pairs.keys():  # Opening bracket
            stack.append(char)
        elif char in pairs.values():  # Closing bracket
            if not stack:
                print("Invalid sequence")
                return
            if char != pairs[stack.pop()]:
                print("Invalid sequence")
                return
    
    # Generate missing closing brackets
    missing = ''
    while stack:
        missing = pairs[stack.pop()] + missing
    
    # Print original + missing
    print(f"Original: {''.join(input_seq)}")
    print(f"Missing: {missing}")
    print(f"Complete sequence: {''.join(input_seq) + missing}")

# Test with the given sequence
sequence = "[ ] < { < { [ ( ) ] } > }"
find_missing_brackets(sequence)