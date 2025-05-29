def complete_parentheses(sequence):
    # Remove the '<' character from the sequence
    sequence = sequence.replace('<', '').strip()
    
    # Initialize a stack to track unmatched opening parentheses
    stack = []
    
    # Iterate through the sequence
    for char in sequence:
        if char == '(':
            stack.append(char)
        elif char == ')':
            if stack:
                stack.pop()
    
    # Calculate the number of closing parentheses needed
    needed_closing = len(stack)
    
    # Construct the completed sequence
    completed_sequence = sequence + ')' * needed_closing
    
    # Format the output
    return f"<<<{completed_sequence}>>>"

# Input sequence
input_sequence = "< ( ( )"

# Get the completed sequence
result = complete_parentheses(input_sequence)

# Print the result
print(result)