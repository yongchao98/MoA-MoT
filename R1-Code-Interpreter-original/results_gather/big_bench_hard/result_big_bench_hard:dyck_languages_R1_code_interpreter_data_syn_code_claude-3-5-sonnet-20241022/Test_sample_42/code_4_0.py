def verify_sequence(opening):
    # Define matching pairs
    pairs = {'{': '}', '<': '>'}
    
    # Stack to store opening brackets
    stack = []
    
    # Process opening sequence
    for char in opening:
        if char in pairs:
            stack.append(char)
    
    # Generate closing sequence in correct order
    result = []
    while stack:
        current = stack.pop()  # get last opening bracket
        result.append(pairs[current])  # add corresponding closing bracket
    
    # Join with spaces to match input format
    final_result = ' '.join(result)
    print(f"Opening sequence: {opening}")
    print(f"Required closing sequence: {final_result}")
    return final_result

# Test with our input
opening_sequence = "{ < <"
closing_sequence = verify_sequence(opening_sequence)