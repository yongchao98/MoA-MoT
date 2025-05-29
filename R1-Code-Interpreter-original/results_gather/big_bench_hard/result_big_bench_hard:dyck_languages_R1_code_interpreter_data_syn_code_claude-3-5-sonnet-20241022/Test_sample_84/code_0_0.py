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

# Original incomplete sequence
incomplete = "( { } ) < { < { } >"

# Complete sequence by adding missing brackets
complete = incomplete + " } >"

print(f"Original sequence: {incomplete}")
print(f"Completed sequence: {complete}")
print(f"Is valid: {is_valid_sequence(complete.replace(' ', ''))}")