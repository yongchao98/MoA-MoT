# Initial sequence
sequence = "[ [ [ [ < { ( ) } ( ) > { ( ) } ] { ( ) }"

# Stack to keep track of opened brackets
stack = []

# Corresponding closing brackets
closing_brackets = {
    '[': ']',
    '{': '}',
    '(': ')',
    '<': '>'
}

# Iterate over the sequence to find the opened brackets
for char in sequence:
    if char in closing_brackets:
        stack.append(char)
    elif char in closing_brackets.values():
        if stack and closing_brackets[stack[-1]] == char:
            stack.pop()

# Generate the closing sequence
closing_sequence = ''.join(closing_brackets[char] for char in reversed(stack))

# Print the final sequence
print(sequence + closing_sequence)