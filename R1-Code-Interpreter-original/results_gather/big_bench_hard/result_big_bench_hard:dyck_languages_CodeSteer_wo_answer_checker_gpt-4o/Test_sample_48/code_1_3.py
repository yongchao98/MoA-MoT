def is_sequence_complete(sequence):
    stack = []
    bracket_map = {')': '(', ']': '[', '}': '{', '>': '<'}
    opening_brackets = set(bracket_map.values())

    for char in sequence:
        if char in opening_brackets:
            stack.append(char)
            print(f"Added '{char}' to stack: {stack}")
        elif char in bracket_map:
            if stack and stack[-1] == bracket_map[char]:
                stack.pop()
                print(f"Matched '{char}', stack after pop: {stack}")
            else:
                print(f"Unmatched '{char}', stack: {stack}")
                return False
        # Ignore spaces and other non-bracket characters
        elif char not in bracket_map.values() and char not in bracket_map.keys():
            continue

    print(f"Final stack state: {stack}")
    return len(stack) == 0

sequence = "[ ] < { < { [ ( ) ] } > }"
if is_sequence_complete(sequence):
    print(f"<<<{sequence}>>>")
else:
    print("The sequence is not complete.")