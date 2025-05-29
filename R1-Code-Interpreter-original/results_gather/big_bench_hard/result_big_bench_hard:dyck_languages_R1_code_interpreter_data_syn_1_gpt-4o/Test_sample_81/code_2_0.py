def is_balanced(sequence):
    stack = []
    matching_bracket = {')': '(', '>': '<', ']': '[', '}': '{'}
    
    for char in sequence:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack or stack[-1] != matching_bracket[char]:
                return False
            stack.pop()
    
    return not stack

# Remove spaces from the sequence
sequence = "(<<>([<{{}}()>])[()][][{}]>"
print(is_balanced(sequence))