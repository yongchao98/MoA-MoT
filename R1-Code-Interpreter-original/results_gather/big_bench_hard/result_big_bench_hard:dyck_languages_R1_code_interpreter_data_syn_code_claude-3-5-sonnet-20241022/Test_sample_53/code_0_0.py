def complete_sequence(sequence):
    stack = []
    pairs = {'{': '}', '(': ')', '<': '>'}
    complete_seq = sequence
    
    # Process existing sequence
    for char in sequence:
        if char in '{(<':
            stack.append(char)
        elif char in '})>':
            if stack and pairs[stack[-1]] == char:
                stack.pop()
    
    # Add closing brackets in reverse order
    while stack:
        complete_seq += pairs[stack.pop()]
    
    print(f"Original sequence: {sequence}")
    print(f"Completed sequence: {complete_seq}")
    return complete_seq

sequence = "{ { } ( ( < > )"
complete_sequence(sequence)