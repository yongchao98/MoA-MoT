def find_missing_brackets(sequence):
    # Given sequence
    input_seq = "[ { ( { } ) } < < ( ) { { < [ { [ ( ) ] } ] > } } > >"
    
    # Remove spaces
    input_seq = input_seq.replace(" ", "")
    
    # Initialize stack
    stack = []
    
    # Define bracket pairs
    brackets = {
        '[': ']',
        '{': '}',
        '(': ')',
        '<': '>'
    }
    
    # Process the sequence
    for char in input_seq:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if not stack:  # Stack is empty but we found closing bracket
                print("Invalid sequence")
                return
            if char == brackets[stack[-1]]:  # Matches with top of stack
                stack.pop()
            else:
                print("Invalid sequence")
                return
    
    # Generate missing closing brackets
    missing = ""
    while stack:
        missing = brackets[stack.pop()] + missing
    
    print(f"Missing brackets: {missing}")
    print(f"Complete sequence: {input_seq + missing}")

find_missing_brackets("[ { ( { } ) } < < ( ) { { < [ { [ ( ) ] } ] > } } > >")