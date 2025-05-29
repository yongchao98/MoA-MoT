def complete_sequence(seq):
    stack = []
    pairs = {'[': ']', '<': '>', '{': '}'}
    for char in seq:
        if char in pairs:
            stack.append(char)
        elif char in pairs.values():
            if not stack or pairs[stack.pop()] != char:
                return "Invalid sequence"
    
    # Add remaining closing brackets
    completed_seq = seq + ''.join(pairs[char] for char in reversed(stack))
    return completed_seq

sequence = "[ [ < < { } >"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)