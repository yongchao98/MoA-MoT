def is_balanced(sequence):
    stack = []
    opening = "({[<"
    closing = ")}]>"
    match = {')': '(', '}': '{', ']': '[', '>': '<'}
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if stack and stack[-1] == match[char]:
                stack.pop()
            else:
                return False
    return not stack

sequence = "( < { < { } > { { ( [ ( { } ) ( ( [ ( [ ( { [ { } ] } ) ] ) ] ) ) [ < [ [ [ [ [ < > ] [ { [ [ { ( ( < [ ] > ) [ ( [ ] ) ] < { [ ] } > ) } ] ] } ] { < ( < > ) > } ] ] ] ] > ] ] ) } } } >"
print(is_balanced(sequence))