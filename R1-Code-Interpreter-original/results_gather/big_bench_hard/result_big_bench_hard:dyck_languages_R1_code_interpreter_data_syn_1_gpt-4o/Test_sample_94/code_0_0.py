def complete_sequence(input_sequence):
    # Define matching pairs
    matching_pairs = {')': '(', ']': '[', '}': '{', '>': '<'}
    opening_brackets = set(matching_pairs.values())
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Iterate through each character in the input sequence
    for char in input_sequence:
        if char in opening_brackets:
            stack.append(char)
        elif char in matching_pairs:
            if stack and stack[-1] == matching_pairs[char]:
                stack.pop()
            else:
                return "The sequence is not properly nested."
    
    # If stack is not empty, it means there are unmatched opening brackets
    if stack:
        return "The sequence is not properly nested."
    
    return "The sequence is properly nested and closed."

# Input sequence
input_sequence = "( ( { < > ( ) [ { < ( < > < ( ( ) ( [ { } ] ) ) > { { } { ( ) } } ) [ ] > } ] { { } } }"

# Check the sequence
result = complete_sequence(input_sequence)
print(result)