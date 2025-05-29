def find_completion(s):
    # Define matching pairs
    pairs = {
        '[': ']',
        '{': '}',
        '(': ')',
        '<': '>'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process the input string
    input_str = s.replace(" ", "")
    
    # First verify if current sequence is valid
    for char in input_str:
        if char in pairs.keys():  # Opening bracket
            stack.append(char)
        elif char in pairs.values():  # Closing bracket
            if not stack:
                return "Invalid sequence"
            if char != pairs[stack.pop()]:
                return "Invalid sequence"
    
    # Generate completion
    completion = ""
    while stack:
        completion = pairs[stack.pop()] + completion
    
    print(f"Original: {input_str}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {input_str + completion}")