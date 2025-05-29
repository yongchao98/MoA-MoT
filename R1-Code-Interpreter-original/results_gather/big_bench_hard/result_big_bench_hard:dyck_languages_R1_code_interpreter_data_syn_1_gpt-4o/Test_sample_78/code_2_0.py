def is_balanced(sequence):
    stack = []
    opening = {'<': '>', '[': ']', '{': '}', '(': ')'}
    closing = {'>': '<', ']': '[', '}': '{', ')': '('}
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if not stack or stack[-1] != closing[char]:
                return False
            stack.pop()
    
    return not stack

sequence = "<[{[{[{(( ))}] }(<<(<>)>{}[[[(<{}>)] ] ]{}>)( {} ){}]"
print(is_balanced(sequence))