# Define the input sequence
input_sequence = "< [ < > ]"

# Initialize a stack to track unmatched opening brackets
stack = []

# Define matching pairs for brackets
matching_brackets = {'<': '>', '[': ']', '{': '}', '(': ')' }

# Iterate over each character in the input sequence
for char in input_sequence:
    if char in matching_brackets:
        # If it's an opening bracket, push to stack
        stack.append(char)
    elif char in matching_brackets.values():
        # If it's a closing bracket, pop from stack if it matches
        if stack and matching_brackets[stack[-1]] == char:
            stack.pop()

# Add the necessary closing brackets to complete the sequence
completed_sequence = input_sequence + ''.join(matching_brackets[char] for char in reversed(stack))

# Output the final answer in the required format
print(f"<<<{completed_sequence}>>>")