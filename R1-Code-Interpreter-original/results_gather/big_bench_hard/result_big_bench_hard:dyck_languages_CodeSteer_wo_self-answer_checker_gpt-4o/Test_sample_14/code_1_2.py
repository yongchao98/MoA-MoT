def complete_sequence(sequence):
    stack = []
    pairs = {'(': ')', '{': '}', '[': ']', '<': '>'}
    
    # Process the sequence
    for char in sequence:
        if char in pairs:  # If it's an opening bracket
            stack.append(char)
        elif char in pairs.values():  # If it's a closing bracket
            if stack and pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Unbalanced sequence"
    
    # Complete the sequence by adding necessary closing brackets
    completed_sequence = sequence
    while stack:
        completed_sequence += pairs[stack.pop()]
    
    return f"<<<{completed_sequence}>>>"

input_sequence = "[ ] { ( ( < { ( ( ( { < > ( ) } ) ) [ { } ] { { ( ) } } ) } < > >"
output = complete_sequence(input_sequence)
print(output)