def complete_sequence(sequence):
    stack = []
    pairs = {'(': ')', '{': '}', '[': ']', '<': '>'}
    for char in sequence:
        if char in pairs:
            stack.append(char)
        elif char in pairs.values():
            if stack and pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Unbalanced sequence"
    
    # If stack is not empty, it means there are unclosed brackets
    while stack:
        sequence += pairs[stack.pop()]
    
    return f"<<<{sequence}>>>"

input_sequence = "[ ] { ( ( < { ( ( ( { < > ( ) } ) ) [ { } ] { { ( ) } } ) } < > >"
output = complete_sequence(input_sequence)
print(output)