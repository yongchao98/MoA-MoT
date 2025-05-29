def complete_sequence(s):
    # Define matching pairs
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # First verify existing sequence is valid
    for char in s:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                return "Invalid sequence"
            if char != pairs[stack.pop()]:
                return "Invalid sequence"
    
    # Add closing brackets for any remaining opening brackets
    completion = ''
    while stack:
        completion += pairs[stack.pop()]
    
    # Print original + completion
    print(f"Original: {s}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {s + completion}")

# Test the input
input_seq = "( [ [ [ { } ] ] { < [ < [ { } ] > ] > }"
complete_sequence(input_seq)