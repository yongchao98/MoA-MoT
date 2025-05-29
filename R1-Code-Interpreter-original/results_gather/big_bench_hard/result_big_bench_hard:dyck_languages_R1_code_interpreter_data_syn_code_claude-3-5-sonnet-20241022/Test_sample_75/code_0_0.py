def find_closing_sequence(input_str):
    # Define matching pairs
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Initialize stack to keep track of opening brackets
    stack = []
    
    # Process input string
    for char in input_str:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                return "Invalid sequence"
            if char != pairs[stack.pop()]:
                return "Invalid sequence"
    
    # Generate closing sequence
    closing_sequence = ''
    while stack:
        closing_sequence = pairs[stack.pop()] + closing_sequence
    
    # Print both input and complete sequence
    complete_sequence = input_str + closing_sequence
    print(f"Complete sequence: {complete_sequence}")
    print(f"Closing sequence needed: {closing_sequence}")

# Input string
input_str = "( < [ < ( ) > [ { [ { ( ) } ] } ] ] > [ { } ] ( < { < [ { } ] > } > ) ) ( ( < ( ) > [ ] [ < { ( ) } > ] [ { } ] ["

find_closing_sequence(input_str)