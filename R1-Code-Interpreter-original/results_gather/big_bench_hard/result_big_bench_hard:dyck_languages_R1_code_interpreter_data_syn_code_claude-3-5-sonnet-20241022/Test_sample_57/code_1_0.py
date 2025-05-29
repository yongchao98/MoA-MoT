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
    
    final_sequence = "< > { }"
    print(f"Final sequence: {final_sequence}")
    print(f"Is balanced: {len(stack) == 0}")
    return len(stack) == 0

verify_sequence("< > { }")