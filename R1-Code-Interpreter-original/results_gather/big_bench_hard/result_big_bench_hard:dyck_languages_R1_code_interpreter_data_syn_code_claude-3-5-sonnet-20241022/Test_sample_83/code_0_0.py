def is_valid_sequence(sequence):
    stack = []
    brackets = {')': '(', '}': '{', ']': '[', '>': '<'}
    
    for char in sequence:
        if char in '({[<':
            stack.append(char)
        elif char in ')}]>':
            if not stack or stack.pop() != brackets[char]:
                return False
    
    # Check if sequence is complete
    if not stack:
        return True
    
    # If incomplete, return the closing sequence needed
    closing_sequence = ''
    while stack:
        bracket = stack.pop()
        if bracket == '{':
            closing_sequence += '}'
        elif bracket == '<':
            closing_sequence += '>'
        elif bracket == '[':
            closing_sequence += ']'
        elif bracket == '(':
            closing_sequence += ')'
    
    print(f"Original sequence: { < [ ]")
    print(f"Needed closing sequence: {closing_sequence}")
    return closing_sequence

input_sequence = "{ < [ ]"
result = is_valid_sequence(input_sequence)