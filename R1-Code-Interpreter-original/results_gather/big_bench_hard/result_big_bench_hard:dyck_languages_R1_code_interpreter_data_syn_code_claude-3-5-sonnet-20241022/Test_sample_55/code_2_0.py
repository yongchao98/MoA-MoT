def verify_sequence(original, completion):
    sequence = original + completion
    stack = []
    pairs = {'(': ')', '<': '>'}
    steps = []
    
    print(f"Original: {original}")
    print(f"Complete: {sequence}")
    
    for char in sequence:
        if char in '(<':  # Opening brackets
            stack.append(char)
            steps.append(f"Added {char}: {stack}")
        else:  # Closing brackets
            if not stack:
                return False
            last_open = stack.pop()
            if pairs[last_open] != char:
                return False
            steps.append(f"Matched {last_open} with {char}: {stack}")
    
    is_valid = len(stack) == 0
    print("\nSteps:")
    for step in steps:
        print(step)
    print(f"\nValid sequence: {is_valid}")
    return is_valid

# Test the sequence
original = "(<<<<>"
completion = ">)>"
verify_sequence(original, completion)