def complete_sequence(s):
    # Stack to keep track of opening brackets
    stack = []
    # Input sequence
    input_seq = s.strip()
    # Mapping of closing to opening brackets
    closing_to_opening = {']': '[', '}': '{', '>': '<', ')': '('}
    opening_to_closing = {'[': ']', '{': '}', '<': '>', '(': ')'}
    
    # Process input sequence
    for char in input_seq:
        if char in '[{<(':  # Opening brackets
            stack.append(char)
        elif char in ']}>)':  # Closing brackets
            if not stack or stack[-1] != closing_to_opening[char]:
                return "Invalid sequence"
            stack.pop()
    
    # Generate completion
    completion = ''
    while stack:
        completion += opening_to_closing[stack.pop()]
    
    # Print original + completion
    print(f"Original: {input_seq}")
    print(f"Completion: {completion}")
    print(f"Complete sequence: {input_seq + completion}")

# Test with the given input
input_seq = "[[<<{[]([]]()])()[]{[((())))<<<>[()]>[]]>>[<>]}>>"
complete_sequence(input_seq)