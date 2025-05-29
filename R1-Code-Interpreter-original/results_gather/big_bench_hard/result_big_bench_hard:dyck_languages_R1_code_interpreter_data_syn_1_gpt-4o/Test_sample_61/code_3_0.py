def is_balanced(sequence):
    stack = []
    opening = {'{', '(', '[', '<'}
    closing = {'}': '{', ')': '(', ']': '[', '>': '<'}
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if stack and stack[-1] == closing[char]:
                stack.pop()
            else:
                return False
    return len(stack) == 0

sequence = "{ < { ( { [ ] } ) } > { } ( [ { { { ( [ < > ] ) } } } ] ) { { ( { } ) } }"
print(is_balanced(sequence))