def find_closing_sequence(input_str):
    # Define matching pairs
    pairs = {
        '[': ']',
        '{': '}',
        '(': ')',
        '<': '>'
    }
    
    # Initialize stack to keep track of opening brackets
    stack = []
    
    # Process input string
    for char in input_str:
        if char in '[{(<':  # Opening brackets
            stack.append(char)
        elif char in ']})>':  # Closing brackets
            if not stack:
                return "Invalid sequence"
            if pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate closing sequence
    closing_sequence = ''
    while stack:
        closing_sequence = pairs[stack.pop()] + closing_sequence
    
    return closing_sequence

# Test with the given input
input_str = "[ ] ( [ [ { < { { ( < > [ ] ) } } < > > } ] ] { }"
input_str = input_str.replace(" ", "")  # Remove spaces
result = find_closing_sequence(input_str)
print(result)