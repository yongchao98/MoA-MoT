def complete_sequence(sequence):
    stack = []
    matching_bracket = {'(': ')', '{': '}', '<': '>'}
    
    for char in sequence:
        if char in matching_bracket:
            stack.append(char)
        elif char in matching_bracket.values():
            if stack and matching_bracket[stack[-1]] == char:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Now, we need to add the missing closing brackets
    missing_closures = ''.join(matching_bracket[char] for char in reversed(stack))
    return sequence + missing_closures

sequence = "< { ( { } ) < > }"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)