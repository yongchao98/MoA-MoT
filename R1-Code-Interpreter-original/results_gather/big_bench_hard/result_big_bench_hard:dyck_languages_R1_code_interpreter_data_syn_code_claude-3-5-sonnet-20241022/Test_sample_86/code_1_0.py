def complete_parentheses(s):
    # Remove spaces and create clean input
    s = s.replace(" ", "")
    
    # Count unclosed parentheses
    stack = []
    for char in s:
        if char == '(':
            stack.append(char)
        elif char == ')':
            if stack:  # if there's an opening bracket to match
                stack.pop()
            else:
                print("Invalid sequence - too many closing brackets")
                return None
    
    # The remaining stack length tells us how many closing brackets we need
    needed_close = len(stack)
    
    # Complete the sequence
    result = s + ')' * needed_close
    
    # Verify the completed sequence
    stack = []
    for char in result:
        if char == '(':
            stack.append(char)
        elif char == ')':
            if stack:
                stack.pop()
            else:
                print("Invalid sequence")
                return None
    
    if not stack:  # if stack is empty, sequence is valid
        print(f"Complete sequence: {result}")
        return result
    else:
        print("Cannot create valid sequence")
        return None

# Test the input
input_str = "( ( ( ) )"
complete_parentheses(input_str)