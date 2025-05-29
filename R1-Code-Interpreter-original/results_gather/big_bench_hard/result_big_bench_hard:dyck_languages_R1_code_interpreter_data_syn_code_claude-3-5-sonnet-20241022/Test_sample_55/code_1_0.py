def verify_sequence(sequence):
    def is_balanced(s):
        stack = []
        brackets = {'(': ')', '<': '>'}
        for char in s:
            if char in brackets:  # Opening bracket
                stack.append(char)
            else:  # Closing bracket
                if not stack:  # Stack is empty
                    return False
                if char != brackets[stack.pop()]:
                    return False
        return len(stack) == 0

    original = "(<<<<>"
    completed = sequence
    
    print(f"Original sequence: {original}")
    print(f"Completed sequence: {completed}")
    print(f"Is balanced: {is_balanced(completed)}")
    
    # Show step by step bracket matching
    stack = []
    for i, char in enumerate(completed):
        if char in '(<':
            stack.append(char)
            print(f"Step {i+1}: Added {char} to stack: {stack}")
        else:
            if stack:
                opening = stack.pop()
                print(f"Step {i+1}: Matched {opening} with {char}, stack: {stack}")

# Test the sequence
verify_sequence("(<<<<>>)>")