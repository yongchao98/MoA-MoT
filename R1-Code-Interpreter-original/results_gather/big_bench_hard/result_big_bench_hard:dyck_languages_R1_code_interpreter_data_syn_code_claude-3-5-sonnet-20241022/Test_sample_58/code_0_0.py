def complete_sequence(sequence):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to map opening brackets to closing brackets
    brackets = {
        '<': '>',
        '(': ')',
        '[': ']',
        '{': '}'
    }
    
    # First, verify the given sequence and build remaining stack
    given = sequence.replace(" ", "")
    for char in given:
        if char in brackets:  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            if not stack or brackets[stack[-1]] != char:
                return "Invalid sequence"
            stack.pop()
    
    # Complete the sequence by adding closing brackets
    result = given
    while stack:
        result += brackets[stack.pop()]
    
    print(f"Complete sequence: {result}")
    return result

# Test the sequence
input_seq = "< ( [ [ [ [ ( ) ] ] { } ] ] )"
complete_sequence(input_seq)