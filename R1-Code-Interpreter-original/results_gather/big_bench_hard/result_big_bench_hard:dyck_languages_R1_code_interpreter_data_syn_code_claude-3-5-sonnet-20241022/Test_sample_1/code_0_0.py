def check_brackets(sequence):
    opening = "< [ ["
    closing = "] ] >"
    stack = []
    
    # First add all opening brackets to stack
    for char in opening:
        if char != ' ':
            stack.append(char)
    
    # Then add closing brackets and check if they match
    result = opening
    for char in closing:
        if char != ' ':
            if not stack:
                return "Invalid sequence"
            if (stack[-1] == '<' and char == '>') or \
               (stack[-1] == '[' and char == ']'):
                stack.pop()
                result += char
            else:
                return "Invalid sequence"
    
    print(f"Valid sequence: {result}")
    return result

check_brackets("< [ [")