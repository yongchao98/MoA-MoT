def is_balanced(s):
    stack = []
    brackets = {'[': ']', '<': '>'}
    opening = set(['[', '<'])
    closing = set([']', '>'])
    
    # Test different possible endings
    test_str = s
    missing = ""
    
    for char in test_str:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if not stack:
                return False, ""
            if char != brackets[stack.pop()]:
                return False, ""
    
    # Add missing closing brackets
    while stack:
        missing += brackets[stack.pop()]
    
    print(f"Original string: {s}")
    print(f"Missing brackets: {missing}")
    print(f"Complete string: {s + missing}")
    return True, missing

# Test the input
input_str = "[ < > ] [ [ < > ]"
is_valid, missing = is_balanced(input_str)