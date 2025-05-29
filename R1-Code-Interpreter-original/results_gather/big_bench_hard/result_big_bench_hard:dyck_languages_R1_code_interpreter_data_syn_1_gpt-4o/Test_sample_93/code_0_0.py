def complete_brackets(sequence):
    stack = []
    for char in sequence:
        if char in "[<":
            stack.append(char)
        elif char in "]>":
            if not stack:
                return "Unbalanced"
            top = stack.pop()
            if (char == ']' and top != '[') or (char == '>' and top != '<'):
                return "Unbalanced"
    
    # Add the necessary closing brackets
    closing_sequence = ""
    while stack:
        top = stack.pop()
        if top == '[':
            closing_sequence += ']'
        elif top == '<':
            closing_sequence += '>'
    
    return sequence + closing_sequence

sequence = "[ < > [ [ < [ ] > ]"
completed_sequence = complete_brackets(sequence)
print(completed_sequence)