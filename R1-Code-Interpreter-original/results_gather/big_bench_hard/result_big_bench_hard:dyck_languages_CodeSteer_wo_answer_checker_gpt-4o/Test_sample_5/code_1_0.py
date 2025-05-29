# Given input sequence
input_sequence = "( [ [ [ { } ] ] { < [ < [ { } ] > ] > }"

# Define matching pairs for brackets
matching_brackets = {')': '(', ']': '[', '}': '{', '>': '<'}

# Stack to keep track of opening brackets
stack = []

# Iterate through each character in the input sequence
for char in input_sequence:
    if char in "([{<":
        # Push opening brackets onto the stack
        stack.append(char)
    elif char in ")]}>":
        # Check if the stack is not empty and matches the corresponding opening bracket
        if stack and stack[-1] == matching_brackets[char]:
            stack.pop()

# Construct the closing sequence needed to balance the input
closing_sequence = ""
# Map opening brackets to their corresponding closing brackets
closing_map = {'(': ')', '[': ']', '{': '}', '<': '>'}
while stack:
    opening_bracket = stack.pop()
    closing_sequence += closing_map[opening_bracket]

# Complete the sequence
complete_sequence = input_sequence + closing_sequence

# Output the final answer
print(f"<<<{complete_sequence}>>>")