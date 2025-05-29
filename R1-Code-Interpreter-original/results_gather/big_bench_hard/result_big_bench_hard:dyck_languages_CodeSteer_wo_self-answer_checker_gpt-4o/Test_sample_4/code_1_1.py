def is_valid_sequence(sequence):
    # Define matching pairs
    matching_pairs = {'>': '<', ')': '(', '}': '{', ']': '['}
    # Initialize an empty stack
    stack = []
    
    for char in sequence:
        if char in '<({[':
            # Push opening brackets onto the stack
            stack.append(char)
        elif char in '>)}]':
            # Check if the stack is not empty and the top of the stack matches
            if stack and stack[-1] == matching_pairs[char]:
                stack.pop()
            else:
                return False
    
    # If the stack is empty, all brackets are properly closed
    return len(stack) == 0

sequence = "< ( ( ( < > ) ) ( { { } [ { } ] [ ] < ( ) > } ) )"
if is_valid_sequence(sequence):
    print(f"<<<{sequence}>>>")
else:
    print("The sequence is not properly closed.")