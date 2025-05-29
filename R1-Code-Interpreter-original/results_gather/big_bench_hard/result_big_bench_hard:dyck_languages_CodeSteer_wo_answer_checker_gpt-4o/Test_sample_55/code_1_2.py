def complete_sequence(sequence):
    stack = []
    for char in sequence:
        if char in '(<':
            stack.append(char)
        elif char == ')' and stack and stack[-1] == '(':
            stack.pop()
        elif char == '>' and stack and stack[-1] == '<':
            stack.pop()

    # Append the necessary closing symbols
    completed_sequence = sequence
    while stack:
        opening = stack.pop()
        if opening == '(':
            completed_sequence += ')'
        elif opening == '<':
            completed_sequence += '>'

    # Output the completed sequence in the specified format
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "( < < < >"
complete_sequence(input_sequence)