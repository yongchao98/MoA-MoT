def complete_sequence(sequence):
    stack = []
    for char in sequence:
        if char in "<(":
            stack.append(char)
        elif char == '>':
            if stack and stack[-1] == '<':
                stack.pop()
            else:
                return "Unbalanced"
        elif char == ')':
            if stack and stack[-1] == '(':
                stack.pop()
            else:
                return "Unbalanced"
    
    # Add the necessary closing brackets
    result = sequence
    while stack:
        top = stack.pop()
        if top == '<':
            result += '>'
        elif top == '(':
            result += ')'
    
    return result

sequence = "< < ( < >"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)