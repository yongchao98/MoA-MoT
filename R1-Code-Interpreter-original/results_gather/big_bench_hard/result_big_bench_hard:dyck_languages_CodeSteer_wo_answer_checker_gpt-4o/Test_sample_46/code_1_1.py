# Input sequence with open and some closed brackets
input_sequence = "[ ( < { [ ( ) ] } >"

# Dictionary to map closing brackets to opening brackets
matching_brackets = {')': '(', ']': '[', '}': '{', '>': '<'}

# Stack to keep track of open brackets
stack = []

# Process the input sequence
for char in input_sequence:
    if char in matching_brackets.values():
        # If it's an opening bracket, push it onto the stack
        stack.append(char)
    elif char in matching_brackets:
        # If it's a closing bracket, pop from the stack if it matches
        if stack and stack[-1] == matching_brackets[char]:
            stack.pop()

# Complete the sequence by adding the necessary closing brackets
completed_sequence = input_sequence + ''.join(matching_brackets[open_bracket] for open_bracket in reversed(stack))

# Output the final answer
print(f"<<<{completed_sequence}>>>")