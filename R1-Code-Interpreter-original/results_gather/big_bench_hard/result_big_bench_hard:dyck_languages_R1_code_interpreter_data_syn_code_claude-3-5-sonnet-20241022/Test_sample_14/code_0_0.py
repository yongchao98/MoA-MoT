def complete_brackets(sequence):
    # Map of closing to opening brackets
    brackets_map = {
        ']': '[',
        '}': '{',
        '>': '<',
        ')': '('
    }
    
    # Map of opening to closing brackets
    opening_to_closing = {
        '[': ']',
        '{': '}',
        '<': '>',
        '(': ')'
    }
    
    # Clean and split the input
    chars = [c for c in sequence if c in '[]{}()<>']
    
    # Stack to track opening brackets
    stack = []
    
    # Process existing sequence
    for char in chars:
        if char in '[{(<':  # Opening bracket
            stack.append(char)
        else:  # Closing bracket
            if stack and stack[-1] == brackets_map[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate closing sequence
    completion = ''
    while stack:
        completion = opening_to_closing[stack.pop()] + completion
    
    # Print original + completion
    print(f"Original: {''.join(chars)}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {''.join(chars) + completion}")
    
# Test the sequence
sequence = "[ ] { ( ( < { ( ( ( { < > ( ) } ) ) [ { } ] { { ( ) } } ) } < > >"
complete_brackets(sequence)