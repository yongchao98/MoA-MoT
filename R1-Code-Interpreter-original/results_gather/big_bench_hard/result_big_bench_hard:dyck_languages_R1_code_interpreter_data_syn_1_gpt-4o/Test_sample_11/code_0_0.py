def is_balanced(sequence):
    stack = []
    matching_bracket = {')': '(', '}': '{', ']': '[', '>': '<'}
    
    for char in sequence:
        if char in matching_bracket.values():
            stack.append(char)
        elif char in matching_bracket.keys():
            if stack == [] or matching_bracket[char] != stack.pop():
                return False
    return stack == []

# Original sequence
sequence = "[{<()>}]"  # Adding the missing closing bracket

# Check if the sequence is balanced
is_balanced_sequence = is_balanced(sequence)
print(is_balanced_sequence)