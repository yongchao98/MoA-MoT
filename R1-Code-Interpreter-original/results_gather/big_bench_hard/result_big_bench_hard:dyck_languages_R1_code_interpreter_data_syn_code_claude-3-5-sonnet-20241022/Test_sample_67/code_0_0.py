def complete_brackets(s):
    # Map of closing brackets to opening brackets
    brackets_map = {')': '(', ']': '[', '}': '{', '>': '<'}
    opening_brackets = {'(', '[', '{', '<'}
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    input_sequence = s.replace(" ", "")
    result = list(input_sequence)
    
    # First pass: check existing sequence and build stack
    for char in input_sequence:
        if char in opening_brackets:
            stack.append(char)
        elif char in brackets_map:
            if not stack or stack[-1] != brackets_map[char]:
                return "Invalid sequence"
            stack.pop()
    
    # Second pass: complete the sequence
    closing_map = {'(': ')', '[': ']', '{': '}', '<': '>'}
    completion = ''
    
    # Add closing brackets in reverse order
    while stack:
        opening_bracket = stack.pop()
        completion += closing_map[opening_bracket]
    
    final_sequence = input_sequence + completion
    
    print(f"Original: {input_sequence}")
    print(f"Completion needed: {completion}")
    print(f"Final sequence: {final_sequence}")
    return completion

# Test the input
input_str = "[[<[([<[]>])]>]({}"
complete_brackets(input_str)