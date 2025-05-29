def is_matching(open_char, close_char):
    pairs = {'{': '}', '[': ']', '(': ')', '<': '>'}
    return pairs.get(open_char) == close_char

def is_valid_sequence(sequence):
    stack = []
    opening = '{[(<'
    closing = '}])>'
    
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if not stack:
                return False
            if not is_matching(stack.pop(), char):
                return False
    return len(stack) == 0

def find_completion(sequence):
    stack = []
    opening = '{[(<'
    closing = '}])>'
    pairs = {'{': '}', '[': ']', '(': ')', '<': '>'}
    
    # Process the sequence
    for char in sequence:
        if char in opening:
            stack.append(char)
        elif char in closing:
            if not stack:
                return "Invalid sequence"
            if not is_matching(stack.pop(), char):
                return "Invalid sequence"
    
    # Generate completion
    completion = ''
    while stack:
        completion += pairs[stack.pop()]
    
    return completion

# Input sequence
input_seq = "{ < { ( { [ ] } ) } > { } ( [ { { { ( [ < > ] ) } } } ] ) { { ( { } ) } }"

# Get completion
completion = find_completion(input_seq)
print(f"Required completion: {completion}")
print(f"Complete sequence: {input_seq + completion}")