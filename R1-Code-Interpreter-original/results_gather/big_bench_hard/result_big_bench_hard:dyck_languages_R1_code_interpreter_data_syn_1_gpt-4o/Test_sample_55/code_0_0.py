def complete_sequence(input_sequence):
    # Initialize counters for opening and closing symbols
    open_paren = 0
    open_angle = 0
    close_paren = 0
    close_angle = 0
    
    # Count the symbols in the input sequence
    for char in input_sequence:
        if char == '(':
            open_paren += 1
        elif char == '<':
            open_angle += 1
        elif char == ')':
            close_paren += 1
        elif char == '>':
            close_angle += 1
    
    # Determine the number of closing symbols needed
    needed_close_paren = open_paren - close_paren
    needed_close_angle = open_angle - close_angle
    
    # Generate the completed sequence
    completed_sequence = input_sequence + ')' * needed_close_paren + '>' * needed_close_angle
    
    return completed_sequence

# Input sequence
input_sequence = "( < < < >"

# Get the completed sequence
completed_sequence = complete_sequence(input_sequence)

# Print the completed sequence
print(completed_sequence)