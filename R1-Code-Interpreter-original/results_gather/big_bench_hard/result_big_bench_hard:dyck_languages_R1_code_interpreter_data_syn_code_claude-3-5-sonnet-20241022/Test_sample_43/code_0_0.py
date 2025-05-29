def check_sequence(sequence):
    # Count opening and closing braces
    opening = sequence.count('{')
    closing = sequence.count('}')
    
    # Calculate how many closing braces needed
    needed = opening - closing
    
    # Complete the sequence
    completed = sequence + ('}' * needed)
    
    # Verify if it's valid
    stack = []
    for char in completed:
        if char == '{':
            stack.append(char)
        elif char == '}':
            if not stack:
                return "Invalid sequence"
            stack.pop()
    
    if stack:
        return "Invalid sequence"
    
    print(f"Original sequence: {sequence}")
    print(f"Completed sequence: {completed}")
    print(f"Number of braces needed: {needed}")

# Test the input
check_sequence("{ { { } }")