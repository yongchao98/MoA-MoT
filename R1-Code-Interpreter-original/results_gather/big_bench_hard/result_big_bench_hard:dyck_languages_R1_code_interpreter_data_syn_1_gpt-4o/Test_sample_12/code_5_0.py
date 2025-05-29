def is_properly_closed(sequence):
    stack = []
    opening = {'(': ')', '[': ']', '{': '}', '<': '>'}
    closing = {')': '(', ']': '[', '}': '{', '>': '<'}
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if not stack or stack[-1] != closing[char]:
                return False
            stack.pop()
    
    return len(stack) == 0

sequence = "[<[(())<()>( {{}}[[[<([])()>]]]{{{{{}}}{{}{<[]>}}}})>)]]>"
print(is_properly_closed(sequence))