def is_balanced(sequence):
    stack = []
    brackets = {'(': ')', '<': '>'}
    
    for char in sequence:
        if char in brackets:  # Opening bracket
            stack.append(char)
        else:  # Closing bracket
            if not stack:  # Stack is empty
                return False
            if char != brackets[stack.pop()]:
                return False
    return len(stack) == 0

# Original incomplete sequence
sequence = "(<<<<>"
# Complete sequence by adding "> > )"
complete_sequence = sequence + ">>)"

print(f"Original sequence: {sequence}")
print(f"Completed sequence: {complete_sequence}")
print(f"Is balanced: {is_balanced(complete_sequence)}")