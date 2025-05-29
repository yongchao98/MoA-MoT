def is_properly_closed(sequence):
    stack = []
    matching_bracket = {')': '(', ']': '[', '}': '{', '>': '<'}
    
    for char in sequence:
        if char in matching_bracket.values():
            stack.append(char)
        elif char in matching_bracket.keys():
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return False
    return not stack

sequence = "< [ { { < ( ) > { < { } > ( < ( ) > { < [ ( { { ( < [ ] > ) } } { ( ( [ [ { } [ ] ] ] ) ) } ) ] > } ) } }"
print(is_properly_closed(sequence))