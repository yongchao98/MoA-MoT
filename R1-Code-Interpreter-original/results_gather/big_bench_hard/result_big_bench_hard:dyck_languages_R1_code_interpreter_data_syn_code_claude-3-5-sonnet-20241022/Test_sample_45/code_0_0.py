def is_valid_sequence(sequence):
    stack = []
    brackets = {')': '(', '}': '{', '>': '<'}
    
    for char in sequence:
        if char in '({<':
            stack.append(char)
        elif char in ')}>' and stack:
            if stack[-1] == brackets[char]:
                stack.pop()
            else:
                return False
    return len(stack) == 0

def complete_sequence(incomplete):
    stack = []
    closing_map = {'(': ')', '{': '}', '<': '>'}
    
    # Build stack of opening brackets
    for char in incomplete:
        if char in '({<':
            stack.append(char)
        elif char in ')}>' and stack:
            if stack[-1] == list(closing_map.keys())[list(closing_map.values()).index(char)]:
                stack.pop()
    
    # Add closing brackets in reverse order
    completion = ''
    while stack:
        completion += closing_map[stack.pop()]
    
    complete = incomplete + completion
    print(f"Complete sequence: {complete}")
    print(f"Is valid: {is_valid_sequence(complete)}")

# Test the incomplete sequence
incomplete = "< { ( ( )"
complete_sequence(incomplete)