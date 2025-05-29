def find_completion(s):
    # Stack to keep track of opening brackets
    stack = []
    # Mapping of closing to opening brackets
    pairs = {'>': '<', '}': '{', ')': '('}
    # Mapping of opening to closing brackets
    closing = {'<': '>', '{': '}', '(': ')'}
    
    # Process existing sequence
    for char in s:
        if char in '<{(':  # Opening brackets
            stack.append(char)
        elif char in '>})':  # Closing brackets
            if not stack or stack[-1] != pairs[char]:
                return "Invalid sequence"
            stack.pop()
    
    # Generate completion
    completion = ''
    # Add closing brackets in reverse order
    while stack:
        completion += closing[stack.pop()]
    
    # Print original + completion and verify it's valid
    result = s + completion
    print(f"Original: {s}")
    print(f"Completion: {completion}")
    print(f"Complete sequence: {result}")
    return completion

# Test the input
input_seq = "< { < > } { ("
find_completion(input_seq)