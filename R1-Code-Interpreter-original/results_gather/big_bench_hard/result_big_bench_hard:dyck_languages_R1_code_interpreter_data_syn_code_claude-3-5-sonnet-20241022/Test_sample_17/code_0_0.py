def complete_sequence(s):
    # Stack to keep track of opening brackets
    stack = []
    # Mapping of closing to opening brackets
    brackets = {')': '(', '>': '<'}
    
    # Process existing sequence
    for char in s:
        if char in '<(':  # Opening brackets
            stack.append(char)
        elif char in ')>':  # Closing brackets
            if stack and stack[-1] == brackets[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Generate completion
    completion = ''
    # Process remaining open brackets in reverse order
    while stack:
        current = stack.pop()
        if current == '(':
            completion += ')'
        elif current == '<':
            completion += '>'
    
    print(f"Original sequence: {s}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {s + completion}")

# Test the input
complete_sequence("< ( ( )")