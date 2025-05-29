def is_valid_with_completion(s):
    stack = []
    pairs = {'}': '{', '>': '<', ')': '('}
    opening = '{(<'
    sequence = s
    needed = []
    
    # Process existing sequence
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in pairs:
            if not stack or stack[-1] != pairs[char]:
                return "Invalid sequence"
            stack.pop()
    
    # Generate completion
    completion = ''
    for bracket in reversed(stack):
        if bracket == '{':
            completion += '}'
        elif bracket == '<':
            completion += '>'
        elif bracket == '(':
            completion += ')'
    
    print(f"Original sequence: {sequence}")
    print(f"Needed completion: {completion}")
    print(f"Complete sequence: {sequence + completion}")
    
# Test the sequence
sequence = "{ } < { } < > ("
is_valid_with_completion(sequence)