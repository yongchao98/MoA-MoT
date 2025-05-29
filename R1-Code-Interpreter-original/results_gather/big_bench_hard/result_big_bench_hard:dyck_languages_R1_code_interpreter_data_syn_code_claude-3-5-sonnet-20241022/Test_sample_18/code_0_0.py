def complete_sequence(s):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to map opening to closing brackets
    brackets = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    # Keep track of missing closing brackets
    missing = []
    
    # Process the input sequence
    for char in s:
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                return "Invalid sequence"
            if char != brackets[stack.pop()]:
                return "Invalid sequence"
    
    # Add missing closing brackets
    while stack:
        missing.append(brackets[stack.pop()])
    
    # Print the complete sequence
    print(f"Original sequence: {s}")
    print(f"Missing brackets: {''.join(missing[::-1])}")
    print(f"Complete sequence: {s + ''.join(missing[::-1])}")

# Test the input
input_seq = "( ( ) ) [ { ( < > ) }"
complete_sequence(input_seq)