def complete_sequence(sequence):
    stack = []
    matching_bracket = {']': '[', '>': '<', ')': '(', '}': '{'}
    opening_brackets = set(matching_bracket.values())
    closing_brackets = set(matching_bracket.keys())
    
    for char in sequence:
        if char in opening_brackets:
            stack.append(char)
        elif char in closing_brackets:
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Now, we need to close the remaining open brackets
    closing_sequence = []
    reverse_matching_bracket = {v: k for k, v in matching_bracket.items()}
    while stack:
        opening = stack.pop()
        closing_sequence.append(reverse_matching_bracket[opening])
    
    return ''.join(closing_sequence)

sequence = "[ [ < [ ( [ < [ ] > ] ) ] > ] ( { }"
closing_sequence = complete_sequence(sequence)
print(closing_sequence)