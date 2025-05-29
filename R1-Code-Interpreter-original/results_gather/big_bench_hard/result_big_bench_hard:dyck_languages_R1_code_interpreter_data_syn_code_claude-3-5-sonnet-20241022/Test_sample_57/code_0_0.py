def is_balanced(sequence):
    stack = []
    pairs = {')': '(', '}': '{', '>': '<'}
    
    for char in sequence:
        if char in '({<':
            stack.append(char)
        elif char in ')}>' and stack:
            if stack[-1] == pairs[char]:
                stack.pop()
            else:
                return False
    
    # Check if sequence is balanced
    sequence_to_test = "< > {"
    print(f"Original sequence: {sequence_to_test}")
    print(f"Stack after processing: {stack}")
    print(f"Need to close with: {['}']}")
    return stack

is_balanced("< > {")