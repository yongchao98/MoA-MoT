def complete_sequence(input_sequence):
    # Define matching pairs
    matching_pairs = {'(': ')', '{': '}'}
    
    # Initialize a stack to keep track of open parentheses and braces
    stack = []
    
    # Iterate through each character in the input sequence
    for char in input_sequence:
        if char in matching_pairs:  # If it's an opening parenthesis or brace
            stack.append(char)
        elif char in matching_pairs.values():  # If it's a closing parenthesis or brace
            if stack and matching_pairs[stack[-1]] == char:
                stack.pop()  # Pop the matching opening from the stack
    
    # Complete the sequence by closing any remaining open parentheses or braces
    completed_sequence = input_sequence
    while stack:
        open_char = stack.pop()
        completed_sequence += matching_pairs[open_char]
    
    # Output the final answer
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "( { ( )"

# Call the function to complete the sequence
complete_sequence(input_sequence)