def complete_sequence(sequence):
    stack = []
    pairs = {'}': '{', ']': '[', ')': '(', '>': '<'}
    opening = set(pairs.values())
    closing = set(pairs.keys())
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if stack and stack[-1] == pairs[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Add the necessary closing brackets
    closing_sequence = ""
    while stack:
        opening_bracket = stack.pop()
        for close, open in pairs.items():
            if open == opening_bracket:
                closing_sequence += close
                break
    
    return sequence + closing_sequence

sequence = "{ { [ { < { } > } ( ) ]"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)