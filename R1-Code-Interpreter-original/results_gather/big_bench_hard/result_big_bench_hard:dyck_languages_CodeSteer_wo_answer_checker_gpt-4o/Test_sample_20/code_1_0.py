def complete_brackets(input_sequence):
    # Define matching pairs of brackets
    matching_brackets = {')': '(', '}': '{', ']': '[', '>': '<'}
    opening_brackets = set(matching_brackets.values())
    
    # Initialize an empty stack
    stack = []
    
    # Iterate through each character in the input sequence
    for char in input_sequence:
        if char in opening_brackets:
            # Push opening brackets onto the stack
            stack.append(char)
        elif char in matching_brackets:
            # Check if the stack is not empty and the top matches
            if stack and stack[-1] == matching_brackets[char]:
                stack.pop()
    
    # Generate the sequence of closing brackets needed
    closing_sequence = ''
    while stack:
        opening_bracket = stack.pop()
        # Find the corresponding closing bracket
        for closing, opening in matching_brackets.items():
            if opening == opening_bracket:
                closing_sequence += closing
                break
    
    # Output the final answer in the required format
    return f"<<<{closing_sequence}>>>"

# Input sequence
input_sequence = "{ ( < [ < > ]"

# Get the final answer
final_answer = complete_brackets(input_sequence)

# Print the final answer
print(final_answer)