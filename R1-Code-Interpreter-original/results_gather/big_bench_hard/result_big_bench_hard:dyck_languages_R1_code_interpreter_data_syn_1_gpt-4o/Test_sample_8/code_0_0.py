def complete_sequence(sequence):
    stack = []
    for char in sequence:
        if char in "[<":
            stack.append(char)
        elif char in "]>":
            if not stack:
                return "Invalid sequence"
            top = stack.pop()
            if (char == ']' and top != '[') or (char == '>' and top != '<'):
                return "Invalid sequence"
    
    # Add remaining closing brackets
    while stack:
        top = stack.pop()
        if top == '[':
            sequence += ']'
        elif top == '<':
            sequence += '>'
    
    return sequence

input_sequence = "[ < > ] [ [ < > ]"
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)