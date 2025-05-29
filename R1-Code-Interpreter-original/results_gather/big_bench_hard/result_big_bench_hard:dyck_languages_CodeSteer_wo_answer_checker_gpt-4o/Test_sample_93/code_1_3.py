def complete_sequence(sequence):
    stack = []
    completed_sequence = []

    for char in sequence:
        if char in "<[":
            stack.append(char)
            completed_sequence.append(char)
        elif char == ">":
            completed_sequence.append(char)
            if stack and stack[-1] == "<":
                stack.pop()
        elif char == "]":
            completed_sequence.append(char)
            if stack and stack[-1] == "[":
                stack.pop()

    # Add the necessary closing brackets for unmatched opening brackets
    while stack:
        char = stack.pop()
        if char == "<":
            completed_sequence.append(">")
        elif char == "[":
            completed_sequence.append("]")

    return ''.join(completed_sequence)

# Given input sequence
input_sequence = "< > [ [ < [ ] > ]"
completed_sequence = complete_sequence(input_sequence)
print(f"<<<{completed_sequence}>>>")