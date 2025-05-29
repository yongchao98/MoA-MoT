def is_valid(sequence):
    stack = []
    pairs = {')': '(', '}': '{', ']': '['}
    
    for char in sequence:
        if char in '({[':
            stack.append(char)
        elif char in ')}]':
            if not stack or stack.pop() != pairs[char]:
                return False
    return len(stack) == 0

def complete_sequence(sequence):
    stack = []
    pairs = {'(': ')', '{': '}', '[': ']'}
    
    # Process the input sequence
    for char in sequence:
        if char in '({[':
            stack.append(char)
        elif char in ')}]':
            if stack and stack[-1] == list(pairs.keys())[list(pairs.values()).index(char)]:
                stack.pop()
    
    # Complete the sequence by adding closing brackets
    completion = ''
    while stack:
        completion += pairs[stack.pop()]
    
    # Return the original sequence plus the completion
    result = sequence + completion
    print(f"Completed sequence: {result}")
    print(f"Is valid: {is_valid(result)}")
    return result

# Test the input
sequence = "( { { } }"
complete_sequence(sequence)