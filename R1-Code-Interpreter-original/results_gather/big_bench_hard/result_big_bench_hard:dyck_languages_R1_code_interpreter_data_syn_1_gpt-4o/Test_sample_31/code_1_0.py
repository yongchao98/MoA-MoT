def complete_sequence(sequence):
    stack = []
    pairs = {'(': ')', '<': '>', '[': ']'}
    for char in sequence:
        if char in pairs:
            stack.append(char)
        elif char in pairs.values():
            if not stack or pairs[stack.pop()] != char:
                return "Invalid sequence"
    
    # Add the necessary closing brackets in the correct order
    closing_sequence = ''.join(pairs[char] for char in reversed(stack))
    return sequence + closing_sequence

input_sequence = "( < [ ( )"
completed_sequence = complete_sequence(input_sequence)
print(completed_sequence)