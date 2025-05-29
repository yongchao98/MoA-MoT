def is_valid(sequence):
    stack = []
    pairs = {'<': '>', '(': ')'}
    
    for char in sequence:
        if char in '<(':  # Opening brackets
            stack.append(char)
        else:  # Closing brackets
            if not stack:  # No opening brackets to match
                return False
            if char != pairs[stack.pop()]:  # Mismatch
                return False
    return len(stack) == 0

def complete_sequence(partial):
    # Convert to list for manipulation
    sequence = list(partial)
    stack = []
    
    # Build stack of opening brackets
    for char in partial:
        if char in '<(':
            stack.append(char)
        elif char in '>)':
            if stack:
                stack.pop()
    
    # Add closing brackets in reverse order
    completion = ''
    pairs = {'<': '>', '(': ')'}
    while stack:
        completion += pairs[stack.pop()]
    
    # Print the complete sequence
    print(f"Original: {partial}")
    print(f"Completion needed: {completion}")
    print(f"Full sequence: {partial + completion}")

# Test the sequence
complete_sequence("<<(<>")