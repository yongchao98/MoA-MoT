def complete_sequence(s):
    # Define matching pairs
    pairs = {'{': '}', '[': ']', '<': '>', '(': ')'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process the input string
    for char in s:
        if char in '{[<(':  # Opening brackets
            stack.append(char)
        elif char in '}]>)':  # Closing brackets
            if not stack:
                return "Invalid sequence"
            if char != pairs[stack.pop()]:
                return "Invalid sequence"
    
    # Generate the completion
    completion = ''
    while stack:
        completion = pairs[stack.pop()] + completion
    
    # Print the full sequence
    print(s + completion)

# Test the input
input_seq = "{ [ < [ < { < { } > } > ( ( < { [ ] } > { { } } ) { } ) ] > ]"
complete_sequence(input_seq)