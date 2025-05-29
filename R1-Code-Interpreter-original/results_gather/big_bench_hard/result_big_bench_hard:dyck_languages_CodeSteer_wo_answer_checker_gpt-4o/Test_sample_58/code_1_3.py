def is_balanced(sequence):
    # Define matching pairs
    matching_pairs = {'(': ')', '[': ']', '{': '}', '<': '>'}
    # Stack to keep track of opening symbols
    stack = []

    for char in sequence:
        if char in matching_pairs:
            # Push opening symbols onto the stack
            stack.append(char)
        elif char in matching_pairs.values():
            # Check if the stack is not empty and the top matches the closing symbol
            if stack and matching_pairs[stack[-1]] == char:
                stack.pop()
            else:
                return False

    # Check if the stack is empty (all symbols are balanced)
    return not stack

# Initial sequence
sequence = "< ( [ [ [ [ ( ) ] ] { } ] ] )"

# Check if the sequence is balanced
# We need to ensure that the sequence starts with '<' and ends with '>'
if is_balanced(sequence) and sequence.startswith('<') and sequence.endswith(')'):
    # Add the closing angle bracket to complete the sequence
    completed_sequence = sequence + " >"
else:
    completed_sequence = "The sequence is not balanced."

# Output the final answer in the required format
print(f"<<<{completed_sequence}>>>")