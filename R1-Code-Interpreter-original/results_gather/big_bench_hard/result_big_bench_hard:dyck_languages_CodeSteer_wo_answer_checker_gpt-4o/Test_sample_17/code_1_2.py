def complete_parentheses(sequence):
    # Remove the '<' character from the sequence
    sequence = sequence.replace('<', '').strip()
    
    # Count the number of opening and closing parentheses
    open_count = sequence.count('(')
    close_count = sequence.count(')')
    
    # Calculate the number of closing parentheses needed
    needed_closing = open_count - close_count
    
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