def is_balanced(sequence):
    stack = []
    opening = "([{<"
    closing = ")]}>"
    matches = {')': '(', ']': '[', '}': '{', '>': '<'}
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if not stack or stack[-1] != matches[char]:
                return False
            stack.pop()
    
    return not stack

sequence = "[ ] < { < { [ ( ) ] } > }"
# Remove spaces for processing
sequence = sequence.replace(" ", "")
print(is_balanced(sequence))