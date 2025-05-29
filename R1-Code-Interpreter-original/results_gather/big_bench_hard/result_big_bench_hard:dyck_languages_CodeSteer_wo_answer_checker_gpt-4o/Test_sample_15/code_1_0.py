# Define matching pairs for opening and closing symbols
matching_pairs = {
    '<': '>',
    '{': '}',
    '(': ')'
}

# Given input sequence
input_sequence = "< { < > } { ("

# Initialize a stack to keep track of unmatched opening symbols
stack = []

# Iterate through each character in the input sequence
for char in input_sequence:
    if char in matching_pairs:  # If it's an opening symbol
        stack.append(char)
    elif char in matching_pairs.values():  # If it's a closing symbol
        if stack and matching_pairs[stack[-1]] == char:
            stack.pop()  # Match found, pop from stack

# Generate the closing sequence for unmatched opening symbols
closing_sequence = ''.join(matching_pairs[stack.pop()] for _ in range(len(stack)))

# Complete sequence
complete_sequence = input_sequence + closing_sequence

# Output the final answer
print(f"<<<{complete_sequence}>>>")