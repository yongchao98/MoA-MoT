def complete_bracket_sequence(sequence):
    # Define matching pairs of brackets
    matching_brackets = {'<': '>', '{': '}', '[': ']', '(': ')'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through the sequence
    for char in sequence:
        if char in matching_brackets:
            # If it's an opening bracket, push to stack
            stack.append(char)
        elif char in matching_brackets.values():
            # If it's a closing bracket, pop from stack
            if stack and matching_brackets[stack[-1]] == char:
                stack.pop()
    
    # Construct the closing sequence
    closing_sequence = ''.join(matching_brackets[stack.pop()] for _ in range(len(stack)))
    
    # Combine the original sequence with the closing sequence
    complete_sequence = sequence + closing_sequence
    
    # Output the final answer
    print(f"<<<{complete_sequence}>>>")

# Input sequence
input_sequence = "< { < [ [ ( { } ) ] ]"
complete_bracket_sequence(input_sequence)