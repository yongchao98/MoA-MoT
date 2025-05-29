def complete_sequence(s):
    # Initialize stack to track opening brackets
    stack = []
    # Original input
    original = s
    
    # First pass: validate existing sequence and build stack
    for char in s:
        if char in ['[', '<']:
            stack.append(char)
        elif char in [']', '>']:
            if not stack:
                return "Invalid sequence"
            if (char == ']' and stack[-1] != '[') or (char == '>' and stack[-1] != '<'):
                return "Invalid sequence"
            stack.pop()
    
    # Build completion
    completion = ''
    while stack:
        bracket = stack.pop()
        if bracket == '[':
            completion += ']'
        elif bracket == '<':
            completion += '>'
    
    # Print original and completed sequence
    print(f"Original: {original}")
    print(f"Completion needed: {completion}")
    print(f"Full sequence: {original + completion}")
    
# Test the sequence
complete_sequence("[ < > [ [ < [ ] > ]")