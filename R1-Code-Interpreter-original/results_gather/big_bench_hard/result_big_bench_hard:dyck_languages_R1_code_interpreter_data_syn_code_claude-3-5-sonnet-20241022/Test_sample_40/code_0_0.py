def complete_brackets(s):
    # Initialize stack to track opening brackets
    stack = []
    # Map of closing to opening brackets
    pairs = {'}': '{', '>': '<', ')': '('}
    
    # Process existing sequence
    for char in s:
        if char in '{(<':
            stack.append(char)
        elif char in '})>':
            if stack and stack[-1] == pairs[char]:
                stack.pop()
            else:
                stack.append(char)

    # Generate completion
    completion = ''
    # Convert stack to list for processing from right to left
    remaining = stack[::-1]
    
    # Process each remaining bracket
    for bracket in remaining:
        if bracket in '})>':
            completion += pairs[bracket]
        elif bracket == '{':
            completion += '}'
        elif bracket == '(':
            completion += ')'
        elif bracket == '<':
            completion += '>'
    
    # Print original sequence and completion
    print(f"Original: {s}")
    print(f"Completion needed: {completion}")
    print(f"Complete sequence: {s + completion}")

# Test the input
input_seq = "{ < { } > { ("
complete_brackets(input_seq)