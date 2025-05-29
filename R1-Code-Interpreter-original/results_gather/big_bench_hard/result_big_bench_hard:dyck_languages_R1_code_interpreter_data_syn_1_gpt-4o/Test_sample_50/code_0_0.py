def complete_sequence(input_sequence):
    stack = []
    closing_brackets = { '(': ')', '[': ']', '<': '>' }
    result = []

    for char in input_sequence:
        if char in closing_brackets:
            stack.append(char)
        elif char in closing_brackets.values():
            if stack and closing_brackets[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    while stack:
        result.append(closing_brackets[stack.pop()])

    return ''.join(result)

input_sequence = "(<>([()"
output = complete_sequence(input_sequence)
print(output)