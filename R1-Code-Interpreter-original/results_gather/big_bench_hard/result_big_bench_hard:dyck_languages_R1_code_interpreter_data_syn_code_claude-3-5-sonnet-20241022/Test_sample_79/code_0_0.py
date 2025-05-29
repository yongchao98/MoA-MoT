def find_matching_brackets(sequence):
    # Define matching pairs
    pairs = {
        '<': '>', 
        '[': ']',
        '(': ')',
        '{': '}'
    }
    
    # Initialize stack to track opening brackets
    stack = []
    
    # Process the input sequence
    for char in sequence:
        if char in '<[({':
            stack.append(char)
        elif char in '>])}':
            if not stack:
                return "Invalid sequence"
            if pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate the closing sequence
    closing_sequence = ''
    while stack:
        closing_sequence = pairs[stack.pop()] + closing_sequence
    
    return closing_sequence

# Input sequence
sequence = "< [ ( [ ( ) [ < [ < > ] < > [ ] > ] < < ( ) ( ) > < { } > > [ < > ] ]"
result = find_matching_brackets(sequence)
print(f"Closing brackets needed: {result}")