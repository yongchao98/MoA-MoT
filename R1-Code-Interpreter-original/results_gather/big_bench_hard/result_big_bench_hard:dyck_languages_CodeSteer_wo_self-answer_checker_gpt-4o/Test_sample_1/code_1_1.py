def complete_sequence(input_sequence):
    # Define matching pairs of brackets
    matching_brackets = {'<': '>', '[': ']', '(': ')', '{': '}'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate over each character in the input sequence
    for char in input_sequence:
        if char in matching_brackets:
            # If it's an opening bracket, push onto the stack
            stack.append(char)
    
    # Construct the closing sequence by popping from the stack
    closing_sequence = ""
    while stack:
        opening_bracket = stack.pop()
        closing_sequence += matching_brackets[opening_bracket]
    
    # Combine the input sequence with the closing sequence
    complete_sequence = input_sequence + closing_sequence
    
    # Output the final answer
    print(f"<<<{complete_sequence}>>>")

# Input sequence
input_sequence = "< [ ["

# Call the function with the input sequence
complete_sequence(input_sequence)