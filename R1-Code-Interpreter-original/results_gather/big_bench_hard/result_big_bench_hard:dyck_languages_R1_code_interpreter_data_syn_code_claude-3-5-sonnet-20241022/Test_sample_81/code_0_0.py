def complete_sequence(s):
    # Stack to keep track of opening brackets
    stack = []
    # Original input
    print("Original:", s)
    
    # Process existing sequence
    for char in s:
        if char in '({[<':
            stack.append(char)
        elif char in ')}]>':
            if not stack:
                return "Invalid sequence"
            if (char == ')' and stack[-1] == '(' or
                char == '}' and stack[-1] == '{' or
                char == ']' and stack[-1] == '[' or
                char == '>' and stack[-1] == '<'):
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate closing sequence
    closing = ''
    while stack:
        last = stack.pop()
        if last == '(': closing += ')'
        elif last == '{': closing += '}'
        elif last == '[': closing += ']'
        elif last == '<': closing += '>'
    
    print("Completion:", closing)
    print("Complete sequence:", s + closing)
    return closing

# Test the input
input_seq = "( < < > ( [ < { { } } ( ) > ] ) [ ( ) ] [ ] [ { } ] >"
result = complete_sequence(input_seq)