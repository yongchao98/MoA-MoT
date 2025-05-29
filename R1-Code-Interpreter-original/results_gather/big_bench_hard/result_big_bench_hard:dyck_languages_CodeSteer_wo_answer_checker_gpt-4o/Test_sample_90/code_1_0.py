def complete_sequence(input_sequence):
    # Define matching pairs
    matching_pairs = {'{': '}', '<': '>', '[': ']', '(': ')'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Traverse the input sequence
    for char in input_sequence:
        if char in matching_pairs:  # If it's an opening bracket
            stack.append(char)
        else:  # If it's a closing bracket
            if stack and matching_pairs[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by adding the necessary closing brackets
    completed_sequence = input_sequence
    while stack:
        completed_sequence += matching_pairs[stack.pop()]
    
    # Output the final answer
    print(f"<<<{completed_sequence}>>>")

# Input sequence
input_sequence = "{ < ( [ ( ( [ < { ( { } ) } < < > > > ] ) ) ] ) >"
complete_sequence(input_sequence)