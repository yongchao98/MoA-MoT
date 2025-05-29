# Initial sequence
sequence = "< ( ( ( [ { } ] )"

# Stack to keep track of opening brackets
stack = []

# Dictionary to map closing brackets to opening brackets
matching_brackets = {')': '(', ']': '[', '}': '{', '>': '<'}

# Result list to build the completed sequence
result = []

# Process each character in the input sequence
for char in sequence:
    if char in matching_brackets.values():  # If it's an opening bracket
        stack.append(char)
    elif char in matching_brackets:  # If it's a closing bracket
        if stack and stack[-1] == matching_brackets[char]:
            stack.pop()
    result.append(char)

# Add the necessary closing brackets for any remaining opening brackets
while stack:
    open_bracket = stack.pop()
    for close_bracket, open_bracket_match in matching_brackets.items():
        if open_bracket_match == open_bracket:
            result.append(close_bracket)
            break

# Output the final answer
completed_sequence = ''.join(result)
print(f"<<<{completed_sequence}>>>")