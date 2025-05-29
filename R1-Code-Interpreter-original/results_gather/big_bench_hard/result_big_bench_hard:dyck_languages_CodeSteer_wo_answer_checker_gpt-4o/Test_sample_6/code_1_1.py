def balance_sequence(sequence):
    stack = []
    # Dictionary to match opening and closing symbols
    matching = {')': '(', '}': '{'}

    # Process each character in the sequence
    for char in sequence:
        if char in '({':
            stack.append(char)
        elif char in ')}':
            if stack and stack[-1] == matching[char]:
                stack.pop()
            else:
                # If there's a mismatch, the sequence is unbalanced
                return "Unbalanced sequence"

    # Add necessary closing symbols for unmatched opening symbols
    completed_sequence = sequence
    while stack:
        open_symbol = stack.pop()
        if open_symbol == '(':
            completed_sequence += ')'
        elif open_symbol == '{':
            completed_sequence += '}'

    return completed_sequence

# Initial sequence
sequence = "( { { } }"

# Get the balanced sequence
balanced_sequence = balance_sequence(sequence)

# Output the completed sequence
print(f"<<<{balanced_sequence}>>>")