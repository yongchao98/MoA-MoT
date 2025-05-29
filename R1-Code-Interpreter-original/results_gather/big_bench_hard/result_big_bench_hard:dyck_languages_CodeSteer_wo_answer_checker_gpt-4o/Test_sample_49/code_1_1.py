def complete_sequence(sequence):
    # Define matching pairs of brackets
    matching_brackets = {']': '[', '}': '{', '>': '<'}
    opening_brackets = set(matching_brackets.values())
    
    # Initialize a stack to track opening brackets
    stack = []
    
    # Iterate through each character in the sequence
    for char in sequence:
        if char in opening_brackets:
            # Push opening brackets onto the stack
            stack.append(char)
        elif char in matching_brackets:
            # Check for matching closing brackets
            if stack and stack[-1] == matching_brackets[char]:
                stack.pop()
            else:
                # If there's a mismatch, add the correct opening bracket
                stack.append(matching_brackets[char])
    
    # Add the necessary closing brackets for unmatched opening brackets
    completed_sequence = sequence + ''.join(matching_brackets[opening] for opening in reversed(stack))
    
    return completed_sequence

# Input sequence
input_sequence = "[ < [ ] > ] [ { [ ] }"

# Get the completed sequence
completed_sequence = complete_sequence(input_sequence)

# Output the final answer in the required format
print(f"<<<{completed_sequence}>>>")