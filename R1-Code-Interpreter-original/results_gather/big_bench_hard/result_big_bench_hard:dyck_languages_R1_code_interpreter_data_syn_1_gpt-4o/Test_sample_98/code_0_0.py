def complete_parentheses(sequence):
    stack = []
    matching_bracket = {')': '(', '}': '{'}
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
    
    # Now, we need to add the closing brackets for the remaining open brackets in the stack
    closing_sequence = ''
    while stack:
        open_bracket = stack.pop()
        for close, open in matching_bracket.items():
            if open == open_bracket:
                closing_sequence += close
    
    return closing_sequence

sequence = "( { ( )"
closing_sequence = complete_parentheses(sequence)
print(closing_sequence)