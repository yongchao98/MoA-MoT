def is_properly_closed(sequence):
    stack = []
    opening = {'(': ')', '[': ']', '{': '}', '<': '>'}
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        else:
            if not stack or opening[stack.pop()] != char:
                return False
    
    return len(stack) == 0

sequence = "[<[(())<()>( {{}}[[[<([])()>]]]{{{{{}}}{{}{<[]>}}}})>)]]>"
print(is_properly_closed(sequence))