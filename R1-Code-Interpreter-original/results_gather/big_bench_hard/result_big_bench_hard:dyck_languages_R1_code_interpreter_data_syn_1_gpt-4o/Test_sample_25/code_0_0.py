def complete_sequence(sequence):
    stack = []
    matching_bracket = {')': '(', '}': '{', ']': '[', '>': '<'}
    opening_bracket = {'(', '{', '[', '<'}
    closing_bracket = {')', '}', ']', '>'}
    
    for char in sequence:
        if char in opening_bracket:
            stack.append(char)
        elif char in closing_bracket:
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Add the necessary closing brackets
    closing_sequence = ""
    while stack:
        open_bracket = stack.pop()
        for close, open in matching_bracket.items():
            if open == open_bracket:
                closing_sequence += close
                break
    
    return sequence + closing_sequence

sequence = "[ < { [ ] } >"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)