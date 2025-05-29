def complete_sequence(s):
    # Define matching brackets
    brackets = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Initialize stack to track opening brackets
    stack = []
    
    # Process existing sequence
    for char in s:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack or brackets[stack[-1]] != char:
                return "Invalid sequence"
            stack.pop()
    
    # Generate completion
    completion = ''
    while stack:
        completion += brackets[stack.pop()]
    
    # Combine original and completion
    result = s + completion
    
    print(f"Original sequence: {s}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {result}")

# Test the sequence
sequence = "[[[[<{()}()>{()}}]{()}]"
complete_sequence(sequence)