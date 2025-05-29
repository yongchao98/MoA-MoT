def verify_sequence(sequence):
    stack = []
    pairs = {'(': ')', '<': '>'}
    valid = True
    steps = []
    
    for i, char in enumerate(sequence):
        if char in '(<':
            stack.append(char)
            steps.append(f"Step {i+1}: Push {char} -> Stack: {stack}")
        else:
            if not stack:
                valid = False
                steps.append(f"Step {i+1}: Error - Closing {char} with empty stack")
                break
            last_open = stack.pop()
            if pairs[last_open] != char:
                valid = False
                steps.append(f"Step {i+1}: Error - Mismatched {last_open} with {char}")
                break
            steps.append(f"Step {i+1}: Match {last_open} with {char} -> Stack: {stack}")
    
    if valid and stack:
        valid = False
        steps.append("Error: Unclosed brackets remain: " + "".join(stack))
    
    print("Testing sequence:", sequence)
    print("\nValidation steps:")
    for step in steps:
        print(step)
    print("\nValid:", valid)
    if valid:
        print("This is a correct sequence!")
    
# Test our proposed solution
original = "(<<<<>"
completion = ">>)"
test_sequence = original + completion
verify_sequence(test_sequence)