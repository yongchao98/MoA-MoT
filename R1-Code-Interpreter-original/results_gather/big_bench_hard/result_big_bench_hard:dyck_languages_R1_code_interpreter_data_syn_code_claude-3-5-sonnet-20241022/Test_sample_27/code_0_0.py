def find_completion(sequence):
    # Map of opening to closing brackets
    brackets = {
        '<': '>',
        '{': '}',
        '(': ')',
        '[': ']'
    }
    
    # Stack to keep track of opening brackets
    stack = []
    
    # Process existing sequence
    for char in sequence:
        if char in brackets.keys():  # Opening bracket
            stack.append(char)
        elif char in brackets.values():  # Closing bracket
            # Find corresponding opening bracket
            for opening, closing in brackets.items():
                if closing == char:
                    if stack and stack[-1] == opening:
                        stack.pop()
                    break
    
    # Generate completion
    completion = ''
    while stack:
        opening = stack.pop()
        completion = brackets[opening] + completion
    
    print(f"Original sequence: {sequence}")
    print(f"Needed completion: {completion}")
    print(f"Complete sequence: {sequence + completion}")

# Test with the given sequence
sequence = "< < { ( < ( ) > ) }"
sequence = sequence.replace(" ", "")
find_completion(sequence)