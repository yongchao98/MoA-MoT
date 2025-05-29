def complete_sequence(sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to match opening and closing brackets
    matching_bracket = {')': '(', ']': '[', '}': '{', '>': '<'}
    # List to store the completed sequence
    completed_sequence = list(sequence)
    
    # Iterate over each character in the sequence
    for char in sequence:
        if char in "([{<":
            # If it's an opening bracket, push to stack
            stack.append(char)
        elif char in ")]}>":
            # If it's a closing bracket, pop from stack
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
    
    # Add the necessary closing brackets in reverse order
    while stack:
        open_bracket = stack.pop()
        for close_bracket, open_bracket_match in matching_bracket.items():
            if open_bracket_match == open_bracket:
                completed_sequence.append(close_bracket)
                break
    
    return ''.join(completed_sequence)

sequence = "( [ [ [ { } ] ] { < [ < [ { } ] > ] > }"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)