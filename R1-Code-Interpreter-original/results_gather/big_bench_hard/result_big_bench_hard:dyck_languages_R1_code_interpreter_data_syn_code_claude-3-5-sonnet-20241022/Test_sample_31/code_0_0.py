def verify_brackets(sequence):
    # Define matching pairs
    pairs = {')': '(', ']': '[', '>': '<'}
    
    # Initial sequence
    initial = "(<[())]>)"
    
    # Stack to track opening brackets
    stack = []
    
    # Check each character
    for char in initial:
        if char in '([<':  # Opening brackets
            stack.append(char)
        elif char in ')]>':  # Closing brackets
            if not stack:  # If stack is empty but we have closing bracket
                return False
            if stack[-1] == pairs[char]:  # If matches with last opening bracket
                stack.pop()
            else:
                return False
    
    # Stack should be empty if all brackets are properly closed
    is_valid = len(stack) == 0
    print(f"Sequence '{initial}' is {'valid' if is_valid else 'invalid'}")
    return is_valid

verify_brackets("(<[())]>)")