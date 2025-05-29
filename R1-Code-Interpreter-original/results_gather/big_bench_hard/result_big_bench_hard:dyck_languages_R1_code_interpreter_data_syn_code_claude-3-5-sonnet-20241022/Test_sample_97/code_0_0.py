def complete_sequence(sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to map closing brackets to opening brackets
    brackets = {')': '(', '>': '<', '}': '{', ']': '['}
    
    # First pass: process existing sequence
    for char in sequence:
        if char in '({[<':
            stack.append(char)
        elif char in ')}]>':
            if not stack:  # If stack is empty but we have closing bracket
                continue
            if stack[-1] == brackets[char]:
                stack.pop()
    
    # Generate completion
    completion = ''
    # Reverse the stack to get proper closing order
    while stack:
        opening = stack.pop()
        # Map opening brackets to closing brackets
        if opening == '(': completion += ')'
        elif opening == '<': completion += '>'
        elif opening == '{': completion += '}'
        elif opening == '[': completion += ']'
    
    print(f"Original sequence: {sequence}")
    print(f"Needed completion: {completion}")
    print(f"Complete sequence: {sequence + completion}")

# Test the sequence
complete_sequence("( < >")