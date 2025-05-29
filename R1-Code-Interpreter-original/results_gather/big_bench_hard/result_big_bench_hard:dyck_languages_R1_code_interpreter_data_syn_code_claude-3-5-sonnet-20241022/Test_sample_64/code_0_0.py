def find_missing_brackets(sequence):
    # Define matching pairs
    pairs = {'[': ']', '<': '>', '(': ')'}
    
    # Convert string to list for easier manipulation
    sequence = list(sequence.replace(" ", ""))
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in sequence:
        if char in '[<(':  # Opening brackets
            stack.append(char)
        elif char in ']>)':  # Closing brackets
            if not stack:
                return "Invalid sequence"
            if char == pairs[stack[-1]]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Add missing closing brackets
    result = sequence.copy()
    while stack:
        result.append(pairs[stack.pop()])
    
    print(''.join(result))

# Test the sequence
sequence = "[ < ( < < > [ ] ( ) >"
find_missing_brackets(sequence)