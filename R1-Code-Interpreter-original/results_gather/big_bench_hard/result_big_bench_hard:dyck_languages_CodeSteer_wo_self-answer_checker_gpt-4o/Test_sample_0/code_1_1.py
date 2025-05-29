# Input sequence
sequence = "[ ["

# Initialize an empty stack
stack = []

# Process each character in the sequence
for char in sequence:
    if char == '[':
        # Push opening bracket onto the stack
        stack.append(char)
    elif char == ']' and stack:
        # Pop from the stack for a matched pair
        stack.pop()

# The stack now contains unmatched opening brackets
# Add the necessary number of closing brackets
completed_sequence = sequence + ']' * len(stack)

# Output the final answer
print(f"<<<{completed_sequence}>>>")