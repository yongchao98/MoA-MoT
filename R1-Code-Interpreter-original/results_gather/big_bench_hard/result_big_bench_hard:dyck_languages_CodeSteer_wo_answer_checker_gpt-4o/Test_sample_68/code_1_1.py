def complete_sequence(input_sequence):
    # Define matching pairs
    matching_pairs = {')': '(', ']': '['}
    opening_symbols = set(matching_pairs.values())
    
    # Initialize a stack to track opening symbols
    stack = []
    
    # Process the input sequence
    for char in input_sequence:
        if char in opening_symbols:
            stack.append(char)
        elif char in matching_pairs:
            if stack and stack[-1] == matching_pairs[char]:
                stack.pop()
            else:
                # If there's a mismatch, the input is invalid
                return "Invalid input sequence"
    
    # Construct the completed sequence
    completed_sequence = input_sequence
    while stack:
        open_symbol = stack.pop()
        if open_symbol == '(':
            completed_sequence += ')'
        elif open_symbol == '[':
            completed_sequence += ']'
    
    return completed_sequence

# Input sequence
input_sequence = "( [ ( )"

# Get the completed sequence
completed_sequence = complete_sequence(input_sequence)

# Output the final answer in the specified format
print(f"<<<{completed_sequence}>>>")