def complete_brackets(s):
    # Stack to keep track of opening brackets
    stack = []
    # Dictionary to match opening and closing brackets
    pairs = {'{': '}', '(': ')', '[': ']'}
    
    # Process existing sequence
    for char in s.replace(" ", ""):
        if char in '{([':
            stack.append(char)
        elif char in '})]':
            if not stack:
                return "Invalid sequence"
            if pairs[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Add closing brackets for remaining opening brackets
    completion = ''
    while stack:
        completion += pairs[stack.pop()]
    
    # Print original + completion
    print(f"Original: {s}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {s + completion}")

# Test the input
input_seq = "{ ( [ [ ] ( ) ]"
complete_brackets(input_seq)