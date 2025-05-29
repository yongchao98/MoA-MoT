def complete_bracket_sequence(sequence):
    # Stack to keep track of opening brackets
    stack = []
    
    # Dictionary to map closing brackets to opening brackets
    bracket_map = {')': '(', ']': '[', '}': '{', '>': '<'}
    
    # Result list to build the completed sequence
    result = []
    
    # Iterate through the sequence
    for char in sequence:
        if char in bracket_map.values():
            # If it's an opening bracket, push to stack
            stack.append(char)
            result.append(char)
        elif char in bracket_map:
            # If it's a closing bracket, check the stack
            if stack and stack[-1] == bracket_map[char]:
                stack.pop()
            result.append(char)
    
    # Close any remaining opening brackets in the stack
    while stack:
        open_bracket = stack.pop()
        # Find the corresponding closing bracket
        for close_bracket, open_bracket_match in bracket_map.items():
            if open_bracket_match == open_bracket:
                result.append(close_bracket)
                break
    
    # Join the result list into a string
    completed_sequence = ''.join(result)
    
    # Output the final answer in the specified format
    print(f"<<<{completed_sequence}>>>")

# Given input sequence
input_sequence = "< { < [ [ ( { } ) ] ]"

# Call the function with the input sequence
complete_bracket_sequence(input_sequence)