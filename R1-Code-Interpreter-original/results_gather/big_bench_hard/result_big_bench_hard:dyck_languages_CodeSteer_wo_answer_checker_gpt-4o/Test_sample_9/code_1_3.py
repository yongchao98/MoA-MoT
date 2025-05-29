def complete_sequence(sequence):
    # Define matching pairs
    matching_bracket = {'[': ']', '(': ')', '{': '}', '<': '>'}
    
    # Stack to keep track of opening brackets and their positions
    stack = []
    
    # Traverse the sequence
    for i, char in enumerate(sequence):
        if char in matching_bracket:  # If it's an opening bracket
            stack.append(char)
        elif char in matching_bracket.values():  # If it's a closing bracket
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by adding the corresponding closing brackets
    completion = []
    while stack:
        open_bracket = stack.pop()
        completion.append(matching_bracket[open_bracket])
    
    # Join the completion to form the final answer
    final_answer = ''.join(completion)
    return f"<<<{final_answer}>>>"

# Input sequence
input_sequence = "[ ] ( [ [ { < { { ( < > [ ] ) } } < > > } ] ] { }"
# Remove spaces for processing
input_sequence = input_sequence.replace(" ", "")

# Get the completed sequence
result = complete_sequence(input_sequence)
print(result)