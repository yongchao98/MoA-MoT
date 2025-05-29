def verify_sequence(original, proposed_completion):
    sequence = original + proposed_completion
    stack = []
    pairs = {'(': ')', '<': '>'}
    steps = []
    
    print(f"Original sequence: {original}")
    print(f"Complete sequence: {sequence}")
    
    for i, char in enumerate(sequence):
        if char in '(<':
            stack.append(char)
            steps.append(f"Step {i+1}: Push {char} -> Stack: {stack}")
        else:
            if not stack:
                steps.append(f"Step {i+1}: Error - Closing {char} with empty stack")
                return False
            last_open = stack.pop()
            if pairs[last_open] != char:
                steps.append(f"Step {i+1}: Error - Mismatched {last_open} with {char}")
                return False
            steps.append(f"Step {i+1}: Match {last_open} with {char} -> Stack: {stack}")
    
    print("\nValidation steps:")
    for step in steps:
        print(step)
    print("\nRemaining stack:", stack)
    return len(stack) == 0

# Test the sequence
original = "(<<<<>"
completion = ">>>)"  # Complete with three '>' and one ')'
verify_sequence(original, completion)