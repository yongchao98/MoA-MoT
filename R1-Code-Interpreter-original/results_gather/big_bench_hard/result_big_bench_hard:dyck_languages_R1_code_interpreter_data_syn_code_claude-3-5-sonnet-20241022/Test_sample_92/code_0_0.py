def find_closing_sequence(s):
    # Define matching pairs
    pairs = {'[': ']', '<': '>', '(': ')'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in s:
        if char in '[<(':
            stack.append(char)
        elif char in ']>)':
            if not stack:
                return "Invalid sequence"
            if char == pairs[stack[-1]]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate closing sequence
    closing = ''
    while stack:
        closing = pairs[stack.pop()] + closing
    
    # Print original sequence + required closing brackets
    print(f"Original: {s}")
    print(f"Needed closing brackets: {closing}")
    print(f"Complete sequence: {s + closing}")

# Test with the given input
input_seq = "[[<[]>]]<"
find_closing_sequence(input_seq)