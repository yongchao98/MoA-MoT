def complete_sequence(sequence):
    stack = []
    matching_bracket = {'>': '<', ']': '[', '}': '{', ')': '('}
    opening_brackets = set(matching_bracket.values())
    
    for char in sequence:
        if char in opening_brackets:
            stack.append(char)
        elif char in matching_bracket:
            if stack and stack[-1] == matching_bracket[char]:
                stack.pop()
            else:
                return "Invalid sequence"
    
    # Now, we need to add the missing closing brackets for the remaining opening brackets in the stack
    closing_sequence = ""
    while stack:
        opening = stack.pop()
        for close, open_ in matching_bracket.items():
            if open_ == opening:
                closing_sequence += close
                break
    
    return sequence + closing_sequence

sequence = "[ < < > > < ( < < < [ [ [ [ ] ] ] ] > > > ) [ < < [ ] > > ] { [ < > ] < > [ ( ) ] } { [ ]"
completed_sequence = complete_sequence(sequence)
print(completed_sequence)