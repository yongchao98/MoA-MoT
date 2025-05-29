def verify_and_fix(sequence):
    stack = []
    correct_sequence = ""
    pairs = {')': '(', '}': '{', '>': '<'}
    
    print("Original sequence:", sequence)
    
    # Process the sequence and build correct version
    for char in sequence:
        if char in '({<':
            stack.append(char)
            correct_sequence += char
        elif char in ')}>' and stack:
            # Check if this is the correct closing bracket for the last opening bracket
            if stack[-1] == pairs[char]:
                stack.pop()
                correct_sequence += char
            else:
                # If not, we need to close the most recent opening bracket first
                last_open = stack[-1]
                if last_open == '(':
                    correct_sequence += ')'
                elif last_open == '{':
                    correct_sequence += '}'
                elif last_open == '<':
                    correct_sequence += '>'
                stack.pop()
    
    print("Corrected sequence:", correct_sequence)
    print("Final stack:", stack)

sequence = "(<{(){}}<>)>"
verify_and_fix(sequence)