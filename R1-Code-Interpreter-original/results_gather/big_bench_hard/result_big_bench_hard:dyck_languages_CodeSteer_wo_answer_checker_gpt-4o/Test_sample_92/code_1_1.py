def complete_sequence(sequence):
    stack = []
    opening_brackets = {'[': ']', '<': '>'}
    closing_brackets = {']': '[', '>': '<'}
    
    # Process each character in the sequence
    for char in sequence:
        if char in opening_brackets:
            stack.append(char)
        elif char in closing_brackets:
            if stack and stack[-1] == closing_brackets[char]:
                stack.pop()
            else:
                # Unbalanced sequence
                return "Unbalanced sequence"
    
    # Add necessary closing brackets to balance the sequence
    while stack:
        open_bracket = stack.pop()
        sequence += opening_brackets[open_bracket]
    
    return f"<<<{sequence}>>>"

# Input sequence
input_sequence = "[ [ < [ ] > ] ] <"
# Output the completed sequence
print(complete_sequence(input_sequence))