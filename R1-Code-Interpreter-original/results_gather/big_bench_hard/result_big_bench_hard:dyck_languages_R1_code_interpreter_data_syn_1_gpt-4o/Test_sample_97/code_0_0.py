def complete_sequence(seq):
    # Initialize a stack to keep track of open brackets
    stack = []
    # Dictionary to match open and close brackets
    matching_bracket = {')': '(', '>': '<'}
    
    # Iterate over each character in the sequence
    for char in seq:
        if char in matching_bracket.values():
            # If it's an open bracket, push to stack
            stack.append(char)
        elif char in matching_bracket.keys():
            # If it's a close bracket, check if it matches the last open bracket
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()  # Pop the matching open bracket
            else:
                # If it doesn't match, the sequence is invalid
                return "Invalid sequence"
    
    # Add the necessary closing brackets to balance the sequence
    completed_sequence = seq + ''.join(matching_bracket[open_bracket] for open_bracket in reversed(stack))
    
    return completed_sequence

# Input sequence
input_sequence = "( < >"
# Get the completed sequence
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)