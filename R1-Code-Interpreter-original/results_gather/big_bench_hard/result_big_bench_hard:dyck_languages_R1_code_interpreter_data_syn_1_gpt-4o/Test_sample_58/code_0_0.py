def complete_sequence(sequence):
    stack = []
    pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
    
    for char in sequence:
        if char in pairs:
            stack.append(char)
        elif char in pairs.values():
            if stack and pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Add the necessary closing characters
    while stack:
        sequence += pairs[stack.pop()]
    
    return sequence

input_sequence = "< ( [ [ [ [ ( ) ] ] { } ] ] )"
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)