def is_balanced(sequence):
    stack = []
    matching_bracket = {')': '(', '}': '{', ']': '['}
    
    for char in sequence:
        if char in matching_bracket.values():
            stack.append(char)
        elif char in matching_bracket.keys():
            if stack == [] or matching_bracket[char] != stack.pop():
                return False
    return stack == []

# Initial sequence
initial_sequence = "( ( { }"

# Complete the sequence
completed_sequence = initial_sequence + " ) )"

# Check if the completed sequence is balanced
print(is_balanced(completed_sequence))