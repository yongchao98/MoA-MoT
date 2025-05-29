def validate_sequence(s):
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    stack = []
    sequence = s.replace(" ", "")
    
    for i, char in enumerate(sequence):
        if char in '([{<':
            stack.append(char)
        elif char in ')]}>':
            if not stack:
                print(f"Error at position {i}: No matching opening bracket for {char}")
                return False
            if pairs[stack[-1]] != char:
                print(f"Error at position {i}: Mismatched brackets. Expected {pairs[stack[-1]]}, found {char}")
                return False
            stack.pop()
    
    if stack:
        print(f"Unclosed brackets remaining: {''.join(stack)}")
        return False
        
    print("Sequence is valid!")
    return True

# Test the completed sequence
sequence = "( ( ( [ < [ < [ ( [ { < ( < ( ( ) ) { } [ ] > ) > } ] ) { < > } ] [ < [ { ( < < { { [ < { [ [ < > [ ] ] ] } > ] { { } } } } > > ) } ] > { < > { } } ] > ( [ ] ) ] > ] ) ( ))))"
validate_sequence(sequence)