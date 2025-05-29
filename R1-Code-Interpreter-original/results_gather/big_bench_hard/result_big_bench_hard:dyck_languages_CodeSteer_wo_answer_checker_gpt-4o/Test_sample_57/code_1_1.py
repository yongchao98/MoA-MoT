def complete_sequence(input_sequence):
    # Define matching pairs
    matching_pairs = {'<': '>', '{': '}'}
    
    # Initialize a stack
    stack = []
    
    # Process the input sequence
    for char in input_sequence:
        if char in matching_pairs:  # If it's an opening bracket/brace
            stack.append(char)
        elif char in matching_pairs.values():  # If it's a closing bracket/brace
            if stack and matching_pairs[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by adding necessary closing brackets/braces
    completed_sequence = input_sequence
    while stack:
        completed_sequence += matching_pairs[stack.pop()]
    
    # Output the final completed sequence
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "< > {"

# Call the function
complete_sequence(input_sequence)