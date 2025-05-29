def complete_brackets(sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to match opening and closing brackets
    brackets = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    
    # Process the input sequence
    input_seq = sequence.replace(" ", "")
    
    # First pass: track opening brackets
    for char in input_seq:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                continue
            if brackets[stack[-1]] == char:
                stack.pop()
    
    # Complete the sequence by adding closing brackets
    completion = ''
    while stack:
        completion += brackets[stack.pop()]
    
    # Print the complete sequence
    print(f"Original: {sequence}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {sequence}{completion}")

# Test the input
complete_brackets("( < > ( [ ( )")