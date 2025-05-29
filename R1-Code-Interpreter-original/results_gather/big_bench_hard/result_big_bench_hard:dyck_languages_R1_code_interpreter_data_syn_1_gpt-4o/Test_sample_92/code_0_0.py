def complete_sequence(sequence):
    stack = []
    closing_brackets = {'[': ']', '<': '>', '{': '}', '(': ')'}
    
    for char in sequence:
        if char in closing_brackets:
            stack.append(closing_brackets[char])
        elif char in closing_brackets.values():
            if stack and stack[-1] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    return ''.join(reversed(stack))

sequence = "[ [ < [ ] > ] ] <"
closing_sequence = complete_sequence(sequence)
print(closing_sequence)