def complete_brackets(sequence):
    # Define matching brackets
    brackets = {
        '{': '}',
        '[': ']',
        '(': ')',
        '<': '>'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process the input sequence
    input_seq = sequence.replace(" ", "")
    
    # First validate and process existing sequence
    for char in input_seq:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if not stack:
                return "Invalid sequence"
            if char != brackets[stack.pop()]:
                return "Invalid sequence"
    
    # Generate closing sequence
    completion = ""
    while stack:
        completion += brackets[stack.pop()]
    
    # Print both the input and completed sequence
    print(f"Input sequence: {input_seq}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {input_seq + completion}")
    return completion

# Test with the given input
input_seq = "{ < { { [ ] } } { < [ { { < > } } [ ( ) ( ) ] [ [ [ [ ( { < ( < ( [ ] ) > ) > } ) ] ] ] ] ] ( ) ( [ ] { } ) > } > [ { ( ( ) ) } ]"
complete_brackets(input_seq)