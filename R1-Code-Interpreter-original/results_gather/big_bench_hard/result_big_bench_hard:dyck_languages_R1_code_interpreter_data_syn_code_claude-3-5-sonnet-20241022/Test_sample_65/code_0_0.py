def complete_brackets(s):
    # Define matching pairs
    pairs = {
        '{': '}',
        '[': ']',
        '<': '>',
        '(': ')'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Input string
    input_str = "{ [ < [ { < ( ) > } ] >"
    
    # Process input string
    for char in input_str:
        if char in '{[<(':
            stack.append(char)
        elif char in '}]>)':
            if not stack:
                return "Invalid sequence"
            if char != pairs[stack.pop()]:
                return "Invalid sequence"
    
    # Generate completion
    completion = ''
    while stack:
        completion = pairs[stack.pop()] + completion
    
    # Print result
    print(f"Original: {input_str}")
    print(f"Completion: {completion}")
    print(f"Complete sequence: {input_str + completion}")

complete_brackets("{ [ < [ { < ( ) > } ] >")