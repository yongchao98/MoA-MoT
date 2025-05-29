def complete_sequence(sequence):
    stack = []
    pairs = {'[': ']', '<': '>', '(': ')', '{': '}'}
    for char in sequence:
        if char in pairs:
            stack.append(char)
        elif char in pairs.values():
            if stack and pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    # Add the necessary closing brackets
    closing_sequence = ''.join(pairs[char] for char in reversed(stack))
    return sequence + closing_sequence

sequence = "[ < < > > < ( < < < [ [ [ [ ] ] ] ] > > > ) [ < < [ ] > > ] { [ < > ] < > [ ( ) ] } { [ ]"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)