def complete_sequence(s):
    # Define matching pairs
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Initialize stack to keep track of opening brackets
    stack = []
    
    # Process the input string
    input_str = s.replace(" ", "")
    
    # First verify the input sequence is valid
    for char in input_str:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>' and stack:
            if pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate the completion
    completion = ''
    # Add closing brackets in reverse order
    while stack:
        completion += pairs[stack.pop()]
    
    print(f"Input: {input_str}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {input_str + completion}")

# Test with the given input
input_str = "(([ << > [ < { } > ] { < [ ] > ( ( ( [ [ { < > } ] ] ) ) ) ( < < > > ) } > ]"
complete_sequence(input_str)