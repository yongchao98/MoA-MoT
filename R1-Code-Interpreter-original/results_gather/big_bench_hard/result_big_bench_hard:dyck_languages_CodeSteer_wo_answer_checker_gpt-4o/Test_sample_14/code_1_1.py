def balance_sequence(sequence):
    stack = []
    matching_bracket = {')': '(', '}': '{', ']': '[', '>': '<'}
    opening_bracket = {'(': ')', '{': '}', '[': ']', '<': '>'}
    
    for char in sequence:
        if char in opening_bracket:
            stack.append(char)
        elif char in matching_bracket:
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                # If there's a mismatch, the sequence is invalid
                return "Invalid sequence"
    
    # Add necessary closing brackets for unmatched opening brackets
    while stack:
        sequence += opening_bracket[stack.pop()]

    return sequence

sequence = "[ ] { ( ( < { ( ( ( { < > ( ) } ) ) [ { } ] { { ( ) } } ) } < > >"
balanced_sequence = balance_sequence(sequence)
print(f"<<<{balanced_sequence}>>>")