def is_valid_sequence(sequence):
    stack = []
    brackets = {'>': '<', '}': '{', ')': '('}
    
    for char in sequence:
        if char in '<{(':
            stack.append(char)
        elif char in '>})':
            if not stack or stack.pop() != brackets[char]:
                return False
    return len(stack) == 0

def complete_sequence(partial):
    original = partial
    stack = []
    # Track opening brackets
    for char in partial:
        if char in '<{(':
            stack.append(char)
        elif char in '>})':
            if stack and stack[-1] == {'>' : '<', '}': '{', ')': '('}[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Add closing brackets for remaining opening brackets
    completion = ''
    while stack:
        bracket = stack.pop()
        if bracket == '<':
            completion += '>'
        elif bracket == '{':
            completion += '}'
        elif bracket == '(':
            completion += ')'
    
    result = original + completion
    print(f"Original sequence: {original}")
    print(f"Completed sequence: {result}")
    print(f"Is valid: {is_valid_sequence(result)}")
    return result

# Test the sequence
sequence = "<<{}>"
complete_sequence(sequence)