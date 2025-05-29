def complete_sequence(sequence):
    stack = []
    for char in sequence:
        if char in '<(':
            stack.append(char)
        elif char == '>':
            if stack and stack[-1] == '<':
                stack.pop()
            else:
                stack.append(char)
        elif char == ')':
            if stack and stack[-1] == '(':
                stack.pop()
            else:
                stack.append(char)

    # Construct the completed sequence
    completed_sequence = sequence
    while stack:
        top = stack.pop()
        if top == '<':
            completed_sequence += '>'
        elif top == '(':
            completed_sequence += ')'

    # Output the completed sequence
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "< < ( < >"
complete_sequence(input_sequence)