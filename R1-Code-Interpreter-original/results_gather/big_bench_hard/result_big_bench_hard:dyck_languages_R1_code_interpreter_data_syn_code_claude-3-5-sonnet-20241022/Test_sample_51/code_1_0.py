def verify_sequence(sequence):
    stack = []
    pairs = {')': '(', '}': '{', '>': '<'}
    
    for char in sequence:
        if char in '({<':
            stack.append(char)
        elif char in ')}>' and stack:
            if stack[-1] == pairs[char]:
                stack.pop()
            else:
                return False
    
    return len(stack) == 0

original = "(<{(){}}<>)>"
complete = "(<{(){}}<>)>)"

print(f"Original sequence balanced: {verify_sequence(original)}")
print(f"Complete sequence balanced: {verify_sequence(complete + '>')}")