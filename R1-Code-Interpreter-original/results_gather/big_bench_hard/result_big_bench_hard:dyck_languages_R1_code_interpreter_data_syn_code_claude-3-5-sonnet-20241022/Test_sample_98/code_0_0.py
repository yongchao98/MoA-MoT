def complete_parentheses(sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to map opening to closing brackets
    brackets = {
        '(': ')',
        '{': '}',
        '[': ']'
    }
    
    # Process the existing sequence
    for char in sequence:
        if char in '({[':
            stack.append(char)
        elif char in ')}]':
            if not stack:
                return "Invalid sequence"
            if brackets[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Complete the sequence by adding closing brackets
    completion = ''
    while stack:
        completion += brackets[stack.pop()]
    
    # Return original sequence + completion
    result = sequence + completion
    print(f"Original sequence: {sequence}")
    print(f"Completed sequence: {result}")
    return result

# Test the sequence
sequence = "({()"
complete_parentheses(sequence)