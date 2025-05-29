# Initial sequence
sequence = "< { < > } { ("

# Stack to track opening brackets
stack = []

# Dictionary to map opening brackets to their corresponding closing brackets
bracket_pairs = {'<': '>', '{': '}', '(': ')' }

# Iterate over each character in the sequence
for char in sequence:
    if char in bracket_pairs:
        # If it's an opening bracket, push to stack
        stack.append(char)
    elif char in bracket_pairs.values():
        # If it's a closing bracket, check if it matches the top of the stack
        if stack and bracket_pairs[stack[-1]] == char:
            stack.pop()

# Add the necessary closing brackets in reverse order
completion = ''.join(bracket_pairs[stack.pop()] for _ in range(len(stack)))

# Output the final completed sequence
final_sequence = sequence + completion
print(f"<<<{final_sequence}>>>")