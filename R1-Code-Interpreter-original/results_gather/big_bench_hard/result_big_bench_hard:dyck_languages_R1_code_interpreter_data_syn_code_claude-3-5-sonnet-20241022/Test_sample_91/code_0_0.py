def complete_brackets(sequence):
    # Define matching brackets
    brackets = {
        '{': '}',
        '[': ']',
        '(': ')',
        '<': '>'
    }
    
    # Initialize stack to keep track of opening brackets
    stack = []
    
    # Process the input sequence
    input_seq = sequence.replace(" ", "")
    result = input_seq
    
    # First pass: check existing sequence and build stack
    for char in input_seq:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if not stack:  # No matching opening bracket
                return "Invalid sequence"
            if char == brackets[stack[-1]]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Add closing brackets for remaining opening brackets
    closing_sequence = ''
    while stack:
        closing_sequence += brackets[stack.pop()]
    
    # Print original sequence and completed sequence
    print(f"Original: {input_seq}")
    print(f"Completed: {input_seq + closing_sequence}")
    return closing_sequence

# Test the input
input_seq = "{{({[[{<[]>}]]})[(<{<>}>)]"
completion = complete_brackets(input_seq)