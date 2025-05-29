def complete_sequence(start):
    # Count opening brackets
    count_square = 0
    for char in start:
        if char == '[':
            count_square += 1
    
    # Generate closing sequence
    closing = ']' * count_square
    
    # Combine and print
    result = start + closing
    print(f"Original: {start}")
    print(f"Completed sequence: {result}")
    print(f"Number of brackets added: {count_square}")
    
# Test the sequence
complete_sequence("[ [")