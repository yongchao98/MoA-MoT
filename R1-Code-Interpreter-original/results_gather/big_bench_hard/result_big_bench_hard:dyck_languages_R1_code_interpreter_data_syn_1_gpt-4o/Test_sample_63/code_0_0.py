def check_brackets(sequence):
    stack = []
    matching_bracket = {']': '[', '>': '<', '}': '{', ')': '('}
    opening_brackets = set(matching_bracket.values())
    
    for char in sequence:
        if char in opening_brackets:
            stack.append(char)
        elif char in matching_bracket:
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # If stack is not empty, it means there are unmatched opening brackets
    if stack:
        return f"Missing closing bracket for {stack[-1]}"
    return "All brackets are properly closed"

sequence = ['[', '<', '<', '{', '}', '>', '<', '>', '>', '(', ')']
result = check_brackets(sequence)
print(result)