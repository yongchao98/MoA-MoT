def complete_sequence(input_sequence):
    # Define matching pairs for brackets
    matching_brackets = {')': '(', '}': '{', '>': '<'}
    opening_brackets = set(matching_brackets.values())
    
    # Initialize a stack to keep track of opening brackets
    stack = []
    
    # Iterate through the input sequence
    for char in input_sequence:
        if char in opening_brackets:
            # Push opening brackets onto the stack
            stack.append(char)
        elif char in matching_brackets:
            # Check if the stack is not empty and matches the current closing bracket
            if stack and stack[-1] == matching_brackets[char]:
                stack.pop()
    
    # Construct the completed sequence
    completed_sequence = input_sequence
    while stack:
        # Pop from the stack and append the corresponding closing bracket
        opening_bracket = stack.pop()
        for closing, opening in matching_brackets.items():
            if opening == opening_bracket:
                completed_sequence += closing
                break
    
    return completed_sequence

# Input sequence
input_sequence = "< ( ( { < > } )"

# Get the completed sequence
completed_sequence = complete_sequence(input_sequence)

# Output the final answer
print(f"<<<{completed_sequence}>>>")