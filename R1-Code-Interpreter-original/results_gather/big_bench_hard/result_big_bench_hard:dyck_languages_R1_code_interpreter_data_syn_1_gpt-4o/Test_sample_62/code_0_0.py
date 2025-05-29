def complete_sequence(sequence):
    # Define matching pairs
    matching_bracket = {')': '(', '>': '<', ']': '[', '}': '{'}
    opening_brackets = set(matching_bracket.values())
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Traverse the sequence
    for char in sequence:
        if char in opening_brackets:
            stack.append(char)
        elif char in matching_bracket:
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # If stack is not empty, add the necessary closing brackets
    closing_sequence = ""
    while stack:
        opening = stack.pop()
        for closing, opening_match in matching_bracket.items():
            if opening_match == opening:
                closing_sequence += closing
                break
    
    return sequence + closing_sequence

# Input sequence
input_sequence = "( ( < [ { { } } ] > { } ) ( < > ) [ ( < > ) ] <"
# Complete the sequence
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)