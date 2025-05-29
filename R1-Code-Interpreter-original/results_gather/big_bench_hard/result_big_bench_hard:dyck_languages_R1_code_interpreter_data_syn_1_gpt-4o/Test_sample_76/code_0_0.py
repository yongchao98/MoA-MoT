def complete_sequence(sequence):
    stack = []
    pairs = {'(': ')', '{': '}', '<': '>'}
    completed_sequence = sequence

    for char in sequence:
        if char in pairs.keys():
            stack.append(char)
        elif char in pairs.values():
            if stack and pairs[stack[-1]] == char:
                stack.pop()

    while stack:
        completed_sequence += pairs[stack.pop()]

    return completed_sequence

sequence = "( { } < > { < { ( < > ) } >"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)