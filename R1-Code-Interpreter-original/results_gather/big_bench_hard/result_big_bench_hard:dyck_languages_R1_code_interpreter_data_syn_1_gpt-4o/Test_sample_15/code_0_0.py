def complete_sequence(input_sequence):
    stack = []
    matching_bracket = {'<': '>', '{': '}', '(': ')'}
    closing_sequence = ''
    
    for char in input_sequence:
        if char in matching_bracket:
            stack.append(char)
        elif char in matching_bracket.values():
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
    
    while stack:
        closing_sequence += matching_bracket[stack.pop()]
    
    return input_sequence + closing_sequence

input_sequence = "< { < > } { ("
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)