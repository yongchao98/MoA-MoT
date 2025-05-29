def is_matching(open_bracket, close_bracket):
    pairs = {
        '(': ')',
        '[': ']',
        '{': '}',
        '<': '>'
    }
    return pairs.get(open_bracket) == close_bracket

def complete_sequence(sequence):
    # Remove spaces
    sequence = sequence.replace(" ", "")
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in sequence:
        if char in '({[<':
            stack.append(char)
        elif char in ')}]>':
            if not stack:
                return "Invalid sequence"
            if not is_matching(stack[-1], char):
                return "Invalid sequence"
            stack.pop()
    
    # Complete the sequence by adding closing brackets
    completion = ''
    while stack:
        bracket = stack.pop()
        if bracket == '(':
            completion += ')'
        elif bracket == '[':
            completion += ']'
        elif bracket == '{':
            completion += '}'
        elif bracket == '<':
            completion += '>'
    
    print(f"Original sequence: {sequence}")
    print(f"Needed completion: {completion}")
    print(f"Complete sequence: {sequence + completion}")

# Test the input sequence
input_seq = "{ ( < > ) } ( ( [ ] ) < [ ( [ [ ] ] [ { } ] { } [ < { [ ] } > ] ( ) ) ]"
complete_sequence(input_seq)