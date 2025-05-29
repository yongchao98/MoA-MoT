def is_balanced(sequence):
    stack = []
    matching_bracket = {')': '(', ']': '[', '>': '<', '}': '{'}
    
    for char in sequence:
        if char in matching_bracket.values():
            stack.append(char)
        elif char in matching_bracket.keys():
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return False, char
    
    return len(stack) == 0, stack

sequence = "[ < [ ] ( ) ( ( { { } } ) ) < { < > } > [ ] > ] < ( ) > ( ( ( ) ) ) ( < >"
is_complete, remaining_stack = is_balanced(sequence)

if not is_complete:
    if remaining_stack:
        # Find the missing closing bracket for the last unclosed opening bracket
        last_open = remaining_stack[-1]
        for close, open in {')': '(', ']': '[', '>': '<', '}': '{'}.items():
            if open == last_open:
                print(close)
                break