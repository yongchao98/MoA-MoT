def complete_sequence(input_sequence):
    stack = []
    matching_bracket = {')': '(', '}': '{', '>': '<'}
    opening_bracket = {'(', '{', '<'}
    closing_bracket = {')', '}', '>'}
    
    for char in input_sequence:
        if char in opening_bracket:
            stack.append(char)
        elif char in closing_bracket:
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Now we need to close any remaining open brackets
    closing_sequence = []
    reverse_matching_bracket = {v: k for k, v in matching_bracket.items()}
    
    while stack:
        opening = stack.pop()
        closing_sequence.append(reverse_matching_bracket[opening])
    
    return ''.join(closing_sequence)

input_sequence = "( < { ( ) { } } ( < > ) >"
closing_sequence = complete_sequence(input_sequence)
print(closing_sequence)