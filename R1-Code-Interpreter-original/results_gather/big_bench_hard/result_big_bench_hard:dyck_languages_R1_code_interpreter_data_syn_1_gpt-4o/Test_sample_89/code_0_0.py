def complete_sequence(sequence):
    stack = []
    opening = {'[': ']', '{': '}', '<': '>', '(': ')'}
    closing = {']': '[', '}': '{', '>': '<', ')': '('}
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if stack and stack[-1] == closing[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Add the necessary closing brackets
    completed_sequence = sequence
    while stack:
        completed_sequence += opening[stack.pop()]
    
    return completed_sequence

sequence = "[ { { < ( ) > } }"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)