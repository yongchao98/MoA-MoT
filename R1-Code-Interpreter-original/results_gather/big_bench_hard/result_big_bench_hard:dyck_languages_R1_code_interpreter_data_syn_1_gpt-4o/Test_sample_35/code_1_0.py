def is_balanced(sequence):
    stack = []
    matching_bracket = {')': '(', '}': '{', ']': '[', '>': '<'}
    
    for char in sequence:
        if char in matching_bracket.values():
            stack.append(char)
        elif char in matching_bracket.keys():
            if stack == [] or matching_bracket[char] != stack.pop():
                return False
        print(f"Current stack: {stack}")
    return stack == []

sequence = "[ { ( { } ) } < < ( ) { { < [ { [ ( ) ] } ] > } } > >"
print(is_balanced(sequence))