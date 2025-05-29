def complete_sequence(seq):
    stack = []
    for char in seq:
        if char in "(<":
            stack.append(char)
        elif char in ")>":
            if not stack:
                return "Unbalanced"
            top = stack.pop()
            if (char == ')' and top != '(') or (char == '>' and top != '<'):
                return "Unbalanced"
    
    # Add remaining closing brackets
    while stack:
        top = stack.pop()
        if top == '(':
            seq += ')'
        elif top == '<':
            seq += '>'
    
    return seq

input_sequence = "( ( < >"
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)