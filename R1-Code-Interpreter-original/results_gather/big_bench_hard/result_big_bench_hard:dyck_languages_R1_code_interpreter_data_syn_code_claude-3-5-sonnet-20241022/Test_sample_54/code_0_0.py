def find_matching_sequence(input_str):
    # Define matching pairs
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Initialize stacks
    stack = []
    result = list(input_str)
    
    # Process input string
    for i, char in enumerate(input_str):
        if char in '({[<':
            stack.append(char)
        elif char in ')}]>':
            if not stack:
                return "Invalid sequence"
            if pairs[stack[-1]] != char:
                return "Invalid sequence"
            stack.pop()
    
    # Complete the sequence
    completion = ''
    while stack:
        completion += pairs[stack.pop()]
    
    # Print original and completion
    print(f"Original sequence: {input_str}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {input_str + completion}")
    
    return completion

# Test the input
input_str = "( < { < { } > { { ( [ ( { } ) ( ( [ ( [ ( { [ { } ] } ) ] ) ] ) ) [ < [ [ [ [ [ < > ] [ { [ [ { ( ( < [ ] > ) [ ( [ ] ) ] < { [ ] } > ) } ] ] } ] { < ( < > ) > } ] ] ] ] > ] ] ) } } } >"
find_matching_sequence(input_str)