def complete_sequence(sequence):
    stack = []
    pairs = {'{': '}', '<': '>'}
    for char in sequence:
        if char in pairs:
            stack.append(char)
    
    completed_sequence = sequence
    while stack:
        completed_sequence += pairs[stack.pop()]
    
    return completed_sequence

sequence = "{ < <"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)