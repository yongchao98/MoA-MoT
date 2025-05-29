# Input sequence
input_sequence = "( < [ ("

# Mapping of opening to closing brackets
bracket_pairs = {
    '(': ')',
    '<': '>',
    '[': ']'
}

# Initialize a stack
stack = []

# Process the input sequence
for char in input_sequence:
    if char in bracket_pairs:
        stack.append(char)

# Generate the closing sequence
closing_sequence = ""
while stack:
    opening_bracket = stack.pop()
    closing_sequence += bracket_pairs[opening_bracket] + " "

# Complete sequence
complete_sequence = input_sequence + " " + closing_sequence.strip()

# Output the final answer
print(f"<<<{complete_sequence}>>>")