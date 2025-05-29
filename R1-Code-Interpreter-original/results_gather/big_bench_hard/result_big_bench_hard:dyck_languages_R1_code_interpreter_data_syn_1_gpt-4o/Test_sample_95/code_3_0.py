def complete_sequence(input_sequence):
    # Define matching pairs
    matching_pairs = {
        ')': '(',
        ']': '[',
        '}': '{',
        '>': '<'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through the input sequence
    for index, char in enumerate(input_sequence):
        if char in matching_pairs.values():
            # If it's an opening bracket, push to stack
            stack.append(char)
        elif char in matching_pairs:
            # If it's a closing bracket, check if it matches the top of the stack
            if stack and stack[-1] == matching_pairs[char]:
                stack.pop()
            else:
                # If it doesn't match, the sequence is not balanced
                return f"Imbalance at index {index}: {char}, Stack: {stack}"
    
    # If the stack is empty, the sequence is balanced
    if not stack:
        return "The sequence is balanced and properly closed"
    else:
        return f"The sequence is not balanced, Remaining Stack: {stack}"

# Adjusted input sequence
adjusted_sequence = "( [ < [ < [ ( [ { < ( < ( ( ) ) { } [ ] > ) > } ] ) { < > } ] [ < [ { ( < < { { [ < { [ [ < > [ ] ] ] } > ] { { } } } } > > ) } ] > { < > { } } ] > ( [ ] ) ] > ] ) ( )"

# Check if the adjusted sequence is balanced
result = complete_sequence(adjusted_sequence)
print(result)