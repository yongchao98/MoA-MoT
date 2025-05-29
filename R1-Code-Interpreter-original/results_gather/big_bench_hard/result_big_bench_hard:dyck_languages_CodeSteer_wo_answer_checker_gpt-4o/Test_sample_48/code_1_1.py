def is_sequence_complete(sequence):
    stack = []
    bracket_map = {')': '(', ']': '[', '}': '{', '>': '<'}
    opening_brackets = set(bracket_map.values())

    for char in sequence:
        if char in opening_brackets:
            stack.append(char)
        elif char in bracket_map:
            if stack and stack[-1] == bracket_map[char]:
                stack.pop()
            else:
                return False

    return len(stack) == 0

sequence = "[ ] < { < { [ ( ) ] } > }"
if is_sequence_complete(sequence):
    print(f"<<<{sequence}>>>")
else:
    print("The sequence is not complete.")