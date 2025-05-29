def complete_sequence(seq):
    stack = []
    pairs = {'[': ']', '<': '>', '(': ')'}
    for char in seq:
        if char in pairs:
            stack.append(char)
        elif char in pairs.values():
            if stack and pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Complete the sequence
    completion = ''.join(pairs[char] for char in reversed(stack))
    return seq + completion

input_sequence = "[ < ( < < > [ ] ( ) >"
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)